#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Seurat)
})

root <- "path omitted"
outdir <- file.path(root, "analysis/05_implantation_failure_atlas/reviewer_stats")
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

heca_centroids <- read.csv(file.path(root, "analysis/05_implantation_failure_atlas/heca_reference_non_endo/tables/heca_non_endo_broad_centroids_logexpr.csv"), check.names = FALSE)
selected_genes <- read.csv(file.path(root, "analysis/05_implantation_failure_atlas/heca_reference_non_endo/tables/heca_non_endo_selected_genes.csv"), check.names = FALSE)$gene
receptivity_sig <- read.csv(file.path(root, "analysis/03_rif/round5_maximal/tables/receptivity_signature_genes_round5.csv"), check.names = FALSE)
late_genes <- na.omit(receptivity_sig$late_genes)
early_genes <- na.omit(receptivity_sig$early_genes)
centroid_mat <- as.matrix(heca_centroids[, -1, drop = FALSE])
rownames(centroid_mat) <- heca_centroids$gene

compute_module_score <- function(expr_mat, pos_genes, neg_genes = character()) {
  up <- intersect(pos_genes, rownames(expr_mat))
  down <- intersect(neg_genes, rownames(expr_mat))
  if (length(up) == 0 && length(down) == 0) return(rep(NA_real_, ncol(expr_mat)))
  z <- t(scale(t(expr_mat)))
  z[!is.finite(z)] <- NA_real_
  up_score <- if (length(up) > 0) colMeans(z[up, , drop = FALSE], na.rm = TRUE) else rep(0, ncol(expr_mat))
  down_score <- if (length(down) > 0) colMeans(z[down, , drop = FALSE], na.rm = TRUE) else rep(0, ncol(expr_mat))
  up_score - down_score
}

project_to_timeline <- function(score, calib_model) {
  pred <- predict(calib_model, newdata = data.frame(receptivity_score = score))
  state <- cut(pred, breaks = c(-Inf, 6.5, 8.5, Inf), labels = c("Delayed", "In-phase", "Advanced"))
  data.frame(pred_day = pred, timing_state = state)
}

read_10x_h5_from_tar <- function(outer_tar, member_name) {
  exdir <- tempfile("tenx_h5_")
  dir.create(exdir)
  utils::untar(outer_tar, exdir = exdir, files = member_name)
  h5_path <- file.path(exdir, member_name)
  Seurat::Read10X_h5(h5_path)
}

read_10x_from_nested_tar <- function(outer_tar, member_name, inner_prefix = NULL) {
  exdir <- tempfile("nested10x_")
  dir.create(exdir)
  utils::untar(outer_tar, exdir = exdir, files = member_name)
  inner_path <- file.path(exdir, member_name)
  utils::untar(inner_path, exdir = exdir)
  files <- list.files(exdir, recursive = TRUE, full.names = TRUE)
  if (!is.null(inner_prefix)) {
    files <- files[grepl(inner_prefix, files) | grepl("filtered_feature_bc_matrix", files)]
  }
  features_path <- files[grepl("features.tsv.gz$", files)][1]
  barcodes_path <- files[grepl("barcodes.tsv.gz$", files)][1]
  matrix_path <- files[grepl("matrix.mtx.gz$", files)][1]
  mtx <- readMM(matrix_path)
  features <- read.delim(gzfile(features_path), header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(gzfile(barcodes_path), header = FALSE, stringsAsFactors = FALSE)
  rownames(mtx) <- make.unique(features[[2]])
  colnames(mtx) <- barcodes[[1]]
  mtx
}

read_10x_flat_from_tar <- function(outer_tar, stem) {
  exdir <- tempfile("flat10x_")
  dir.create(exdir)
  members <- c(
    paste0(stem, "_features.tsv.gz"),
    paste0(stem, "_barcodes.tsv.gz"),
    paste0(stem, "_matrix.mtx.gz")
  )
  utils::untar(outer_tar, exdir = exdir, files = members)
  mtx <- readMM(file.path(exdir, paste0(stem, "_matrix.mtx.gz")))
  features <- read.delim(gzfile(file.path(exdir, paste0(stem, "_features.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)
  barcodes <- read.delim(gzfile(file.path(exdir, paste0(stem, "_barcodes.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)
  rownames(mtx) <- make.unique(features[[2]])
  colnames(mtx) <- barcodes[[1]]
  mtx
}

normalize_selected <- function(mtx, genes) {
  keep <- intersect(genes, rownames(mtx))
  expr <- mtx[keep, , drop = FALSE]
  libsize <- Matrix::colSums(mtx)
  norm <- log1p(t(t(expr) / pmax(libsize, 1)) * 1e4)
  as.matrix(norm)
}

assign_cells_heca <- function(mtx, centroid_mat, min_cor = 0.03) {
  common <- intersect(rownames(centroid_mat), rownames(mtx))
  expr <- normalize_selected(mtx, common)
  ref <- centroid_mat[common, , drop = FALSE]
  cors <- suppressWarnings(cor(expr, ref, use = "pairwise.complete.obs"))
  cors[!is.finite(cors)] <- -1
  best_idx <- max.col(cors, ties.method = "first")
  best_ct <- colnames(cors)[best_idx]
  best_cor <- cors[cbind(seq_len(nrow(cors)), best_idx)]
  best_ct[best_cor < min_cor | !is.finite(best_cor)] <- "Unknown"
  data.frame(
    barcode = colnames(expr),
    assigned_celltype = best_ct,
    assignment_cor = best_cor,
    stringsAsFactors = FALSE
  )
}

aggregate_target_counts <- function(mtx, assignments, target_genes) {
  keep <- intersect(target_genes, rownames(mtx))
  expr <- mtx[keep, , drop = FALSE]
  out <- lapply(unique(assignments$assigned_celltype), function(ct) {
    bc <- assignments$barcode[assignments$assigned_celltype == ct]
    common_bc <- intersect(colnames(expr), bc)
    if (length(common_bc) == 0) return(NULL)
    counts <- Matrix::rowSums(expr[, common_bc, drop = FALSE])
    data.frame(
      gene = keep,
      celltype = ct,
      value = as.numeric(counts),
      n_cells = length(common_bc),
      mean_cor = mean(assignments$assignment_cor[match(common_bc, assignments$barcode)], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  bind_rows(out)
}

bootstrap_diff <- function(case, ctrl, nboot = 2000) {
  obs <- mean(case, na.rm = TRUE) - mean(ctrl, na.rm = TRUE)
  boots <- replicate(nboot, {
    mean(sample(case, length(case), replace = TRUE), na.rm = TRUE) -
      mean(sample(ctrl, length(ctrl), replace = TRUE), na.rm = TRUE)
  })
  c(
    obs = obs,
    lower = as.numeric(quantile(boots, 0.025, na.rm = TRUE)),
    upper = as.numeric(quantile(boots, 0.975, na.rm = TRUE))
  )
}

read_series_meta <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  parse_line <- function(prefix) {
    line <- grep(paste0("^", prefix), lines, value = TRUE)[1]
    vals <- strsplit(line, "\t")[[1]]
    vals <- gsub('^"|"$', "", vals)
    vals[-1]
  }
  meta <- data.frame(
    title = parse_line("!Sample_title"),
    geo_accession = parse_line("!Sample_geo_accession"),
    stringsAsFactors = FALSE
  )
  char_lines <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)
  if (length(char_lines) > 0) {
    for (i in seq_along(char_lines)) {
      vals <- strsplit(char_lines[i], "\t")[[1]]
      vals <- gsub('^"|"$', "", vals)
      meta[[paste0("c", i)]] <- vals[-1]
    }
  }
  meta
}

calc_scores_from_agg <- function(agg_df, meta_cols = c("sample", "group", "day")) {
  celltypes <- sort(setdiff(unique(agg_df$celltype), "Unknown"))
  out <- list()
  for (ct in celltypes) {
    sub <- agg_df %>% filter(celltype == ct)
    mat <- sub %>% select(gene, sample, value) %>% pivot_wider(names_from = sample, values_from = value, values_fill = 0)
    mat <- as.data.frame(mat)
    rownames(mat) <- mat$gene
    mat$gene <- NULL
    mat <- log2(as.matrix(mat) + 1)
    tmp <- sub %>% distinct(across(all_of(meta_cols)))
    score_df <- data.frame(
      sample = colnames(mat),
      celltype = ct,
      receptivity_score = compute_module_score(mat, late_genes, early_genes),
      stringsAsFactors = FALSE
    ) %>% left_join(tmp, by = "sample")
    out[[ct]] <- score_df
  }
  bind_rows(out)
}

load_or_build <- function(path, builder) {
  if (file.exists(path)) return(readRDS(path))
  obj <- builder()
  saveRDS(obj, path)
  obj
}

process_rif <- function() {
  meta_rif <- read.csv(file.path(root, "analysis/03_rif/metadata/GSE250130_metadata_round2.csv"), check.names = FALSE)
  meta_rif$day <- dplyr::case_when(
    meta_rif$group_simple == "LH3" ~ 3,
    meta_rif$group_simple == "LH5" ~ 5,
    meta_rif$group_simple == "Fertile_LH7" ~ 7,
    meta_rif$group_simple == "LH9" ~ 9,
    meta_rif$group_simple == "LH11" ~ 11,
    TRUE ~ NA_real_
  )
  members <- utils::untar("path/to/local/external-storage/03_rif/GSE250130/GSE250130_RAW.tar", list = TRUE)
  member_map <- data.frame(geo_accession = sub("_.*$", "", members), member = members, stringsAsFactors = FALSE)
  meta_rif <- left_join(meta_rif, member_map, by = "geo_accession")
  out <- lapply(seq_len(nrow(meta_rif)), function(i) {
    mtx <- read_10x_from_nested_tar("path/to/local/external-storage/03_rif/GSE250130/GSE250130_RAW.tar", meta_rif$member[i])
    assign <- assign_cells_heca(mtx, centroid_mat)
    agg <- aggregate_target_counts(mtx, assign, selected_genes)
    agg$sample <- meta_rif$sample[i]
    agg$group <- meta_rif$group_simple[i]
    agg$day <- meta_rif$day[i]
    assign_sum <- assign %>% summarise(sample = meta_rif$sample[i], group = meta_rif$group_simple[i], mean_cor = mean(assignment_cor, na.rm = TRUE))
    list(agg = agg, sum = assign_sum)
  })
  list(agg = bind_rows(lapply(out, `[[`, "agg")), summary = bind_rows(lapply(out, `[[`, "sum")))
}

process_endo179640 <- function() {
  meta <- read.csv(file.path(root, "analysis/02_endometriosis/metadata/GSE179640_metadata_round2.csv"), check.names = FALSE)
  meta <- meta %>% filter(!grepl("bulk", supplementary_file_1, ignore.case = TRUE))
  members <- utils::untar("path/to/local/external-storage/reference_annotations/GSE179640/raw/GSE179640_RAW.tar", list = TRUE)
  h5_members <- members[grepl("\\.h5$", members)]
  meta$member <- basename(meta$supplementary_file_1)
  meta <- meta %>% filter(member %in% h5_members)
  out <- lapply(seq_len(nrow(meta)), function(i) {
    mtx <- read_10x_h5_from_tar("path/to/local/external-storage/reference_annotations/GSE179640/raw/GSE179640_RAW.tar", meta$member[i])
    assign <- assign_cells_heca(mtx, centroid_mat)
    agg <- aggregate_target_counts(mtx, assign, selected_genes)
    agg$sample <- meta$geo_accession[i]
    agg$group <- meta$analysis_group[i]
    assign_sum <- assign %>% summarise(sample = meta$geo_accession[i], group = meta$analysis_group[i], mean_cor = mean(assignment_cor, na.rm = TRUE))
    list(agg = agg, sum = assign_sum)
  })
  list(agg = bind_rows(lapply(out, `[[`, "agg")), summary = bind_rows(lapply(out, `[[`, "sum")))
}

process_endo214411 <- function() {
  meta_a <- read_series_meta("path/to/local/external-storage/reference_annotations/GSE214411/metadata/GSE214411-GPL24676_series_matrix.txt.gz")
  meta_a$disease <- sub("^disease: ", "", meta_a$c2)
  meta_a$phase <- sub("^cycle phase: ", "", meta_a$c3)
  meta_b <- read_series_meta("path/to/local/external-storage/reference_annotations/GSE214411/metadata/GSE214411-GPL11154_series_matrix.txt.gz")
  meta_b$disease <- "Control"
  meta_b$phase <- NA_character_
  meta <- bind_rows(meta_a, meta_b)
  members <- utils::untar("path/to/local/external-storage/reference_annotations/GSE214411/raw/GSE214411_RAW.tar", list = TRUE)
  stems <- unique(sub("_(features|barcodes|matrix)\\.tsv\\.gz$|_matrix\\.mtx\\.gz$", "", members))
  stem_df <- data.frame(stem = stems, geo_accession = sub("_.*$", "", stems), stringsAsFactors = FALSE)
  meta <- left_join(meta, stem_df, by = "geo_accession")
  out <- lapply(seq_len(nrow(meta)), function(i) {
    mtx <- read_10x_flat_from_tar("path/to/local/external-storage/reference_annotations/GSE214411/raw/GSE214411_RAW.tar", meta$stem[i])
    assign <- assign_cells_heca(mtx, centroid_mat)
    agg <- aggregate_target_counts(mtx, assign, selected_genes)
    agg$sample <- meta$geo_accession[i]
    agg$group <- meta$disease[i]
    agg$phase <- meta$phase[i]
    assign_sum <- assign %>% summarise(sample = meta$geo_accession[i], group = meta$disease[i], phase = meta$phase[i], mean_cor = mean(assignment_cor, na.rm = TRUE))
    list(agg = agg, sum = assign_sum)
  })
  list(agg = bind_rows(lapply(out, `[[`, "agg")), summary = bind_rows(lapply(out, `[[`, "sum")))
}

rif <- load_or_build(file.path(tabdir, "rif_heca_non_endo_raw_cache.rds"), process_rif)
endo179 <- load_or_build(file.path(tabdir, "endo179640_heca_non_endo_raw_cache.rds"), process_endo179640)
endo214 <- load_or_build(file.path(tabdir, "endo214411_heca_non_endo_raw_cache.rds"), process_endo214411)

full_assign <- read.csv(file.path(root, "analysis/05_implantation_failure_atlas/heca_upgrade/tables/heca_assignment_summary.csv"), check.names = FALSE)
non_endo_assign <- bind_rows(
  rif$summary %>% mutate(dataset = "GSE250130"),
  endo179$summary %>% mutate(dataset = "GSE179640"),
  endo214$summary %>% mutate(dataset = "GSE214411")
) %>%
  group_by(dataset) %>%
  summarise(mean_cor_non_endo = mean(mean_cor, na.rm = TRUE), .groups = "drop")
full_assign_sum <- full_assign %>%
  group_by(dataset) %>%
  summarise(mean_cor_full = mean(mean_cor, na.rm = TRUE), .groups = "drop")
assign_compare <- full_assign_sum %>%
  inner_join(non_endo_assign, by = "dataset") %>%
  mutate(delta = mean_cor_non_endo - mean_cor_full)
write.csv(assign_compare, file.path(tabdir, "heca_non_endo_assignment_quality_comparison.csv"), row.names = FALSE)

rif_scores <- calc_scores_from_agg(rif$agg, c("sample", "group", "day"))
fertile_ref <- rif_scores %>% filter(!is.na(day))
celltype_models <- lapply(split(fertile_ref, fertile_ref$celltype), function(df) {
  if (nrow(df) >= 5 && length(unique(df$day)) >= 3) lm(day ~ receptivity_score, data = df) else NULL
})
global_model <- lm(day ~ receptivity_score, data = fertile_ref)
project_external <- function(score_df) {
  bind_rows(lapply(split(score_df, score_df$celltype), function(df) {
    fit <- celltype_models[[unique(df$celltype)]]
    if (is.null(fit)) fit <- global_model
    proj <- project_to_timeline(df$receptivity_score, fit)
    bind_cols(df, proj)
  }))
}
rif_scores <- project_external(rif_scores)
endo179_scores <- project_external(calc_scores_from_agg(endo179$agg, c("sample", "group")))
endo214_scores <- project_external(calc_scores_from_agg(endo214$agg, c("sample", "group", "phase")))

rif_stats <- rif_scores %>%
  filter(group %in% c("Fertile_LH7", "RIF"), celltype == "Glandular_Epi") %>%
  summarise(
    estimate = bootstrap_diff(pred_day[group == "RIF"], pred_day[group == "Fertile_LH7"])[["obs"]],
    lower = bootstrap_diff(pred_day[group == "RIF"], pred_day[group == "Fertile_LH7"])[["lower"]],
    upper = bootstrap_diff(pred_day[group == "RIF"], pred_day[group == "Fertile_LH7"])[["upper"]]
  ) %>%
  mutate(dataset = "GSE250130", comparison = "RIF vs Fertile_LH7", celltype = "Glandular_Epi")

endo179_stats <- endo179_scores %>%
  filter(group %in% c("Control_Eutopic", "Endo_Eutopic"), celltype == "Stroma") %>%
  summarise(
    estimate = bootstrap_diff(pred_day[group == "Endo_Eutopic"], pred_day[group == "Control_Eutopic"])[["obs"]],
    lower = bootstrap_diff(pred_day[group == "Endo_Eutopic"], pred_day[group == "Control_Eutopic"])[["lower"]],
    upper = bootstrap_diff(pred_day[group == "Endo_Eutopic"], pred_day[group == "Control_Eutopic"])[["upper"]]
  ) %>%
  mutate(dataset = "GSE179640", comparison = "Endo_Eutopic vs Control_Eutopic", celltype = "Stroma")

endo214_stats <- endo214_scores %>%
  filter(group %in% c("Control", "endometriosis"), celltype == "Endothelial") %>%
  summarise(
    estimate = bootstrap_diff(pred_day[group == "endometriosis"], pred_day[group == "Control"])[["obs"]],
    lower = bootstrap_diff(pred_day[group == "endometriosis"], pred_day[group == "Control"])[["lower"]],
    upper = bootstrap_diff(pred_day[group == "endometriosis"], pred_day[group == "Control"])[["upper"]]
  ) %>%
  mutate(dataset = "GSE214411", comparison = "Endometriosis vs Control", celltype = "Endothelial")

non_endo_key <- bind_rows(rif_stats, endo179_stats, endo214_stats) %>%
  select(dataset, comparison, celltype, estimate_non_endo = estimate, lower_non_endo = lower, upper_non_endo = upper)
full_key <- read.csv(file.path(root, "analysis/05_implantation_failure_atlas/heca_upgrade/tables/heca_key_celltype_shift_summary.csv"), check.names = FALSE) %>%
  filter(
    (dataset == "GSE250130" & celltype == "Glandular_Epi") |
      (dataset == "GSE179640" & celltype == "Stroma") |
      (dataset == "GSE214411" & celltype == "Endothelial")
  ) %>%
  select(dataset, comparison, celltype, estimate_full = estimate, lower_full = lower, upper_full = upper)

key_compare <- full_key %>%
  inner_join(non_endo_key, by = c("dataset", "comparison", "celltype")) %>%
  mutate(delta = estimate_non_endo - estimate_full)
write.csv(key_compare, file.path(tabdir, "heca_non_endo_key_shift_comparison.csv"), row.names = FALSE)

p_assign <- assign_compare %>%
  pivot_longer(cols = c(mean_cor_full, mean_cor_non_endo), names_to = "reference", values_to = "mean_cor") %>%
  mutate(reference = recode(reference, mean_cor_full = "Full HECA", mean_cor_non_endo = "Non-endometriosis HECA")) %>%
  ggplot(aes(dataset, mean_cor, color = reference, group = dataset)) +
  geom_line(color = "grey75", linewidth = 0.6) +
  geom_point(size = 3) +
  theme_bw(base_size = 12) +
  labs(title = "Assignment quality under full versus non-endometriosis HECA references", x = "", y = "Mean assignment correlation", color = "")
ggsave(file.path(figdir, "Figure_ATLAS_RS_19_heca_non_endo_assignment_compare.png"), p_assign, width = 7.4, height = 5.6, dpi = 260)

plot_df <- bind_rows(
  key_compare %>% transmute(label = paste(dataset, celltype, sep = " | "), reference = "Full HECA", estimate = estimate_full, lower = lower_full, upper = upper_full),
  key_compare %>% transmute(label = paste(dataset, celltype, sep = " | "), reference = "Non-endometriosis HECA", estimate = estimate_non_endo, lower = lower_non_endo, upper = upper_non_endo)
) %>%
  mutate(label = factor(label, levels = rev(unique(paste(key_compare$dataset, key_compare$celltype, sep = " | ")))))

p_shift <- ggplot(plot_df, aes(estimate, label, color = reference)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), position = position_dodge(width = 0.55), height = 0.2) +
  geom_point(position = position_dodge(width = 0.55), size = 3) +
  theme_bw(base_size = 12) +
  labs(title = "Key cell-state shifts under full versus non-endometriosis HECA references", x = "Predicted-day shift", y = "", color = "")
ggsave(file.path(figdir, "Figure_ATLAS_RS_20_heca_non_endo_shift_compare.png"), p_shift, width = 8.2, height = 4.8, dpi = 260)