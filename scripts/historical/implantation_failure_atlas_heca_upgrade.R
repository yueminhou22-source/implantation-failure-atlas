#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
  library(Seurat)
  library(readxl)
})

root <- "path omitted"
outdir <- file.path(root, "analysis/05_implantation_failure_atlas/heca_upgrade")
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

heca_centroids <- read.csv(file.path(root, "analysis/05_implantation_failure_atlas/heca_reference/tables/heca_broad_centroids_logexpr.csv"), check.names = FALSE)
heca_stage <- read.csv(file.path(root, "analysis/05_implantation_failure_atlas/heca_reference/tables/heca_broad_celltype_stage_counts.csv"), check.names = FALSE)
selected_genes <- read.csv(file.path(root, "analysis/05_implantation_failure_atlas/heca_reference/tables/heca_selected_genes.csv"), check.names = FALSE)$gene

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

bootstrap_diff <- function(case, ctrl, nboot = 3000) {
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

perm_p <- function(x, g, nperm = 3000) {
  obs <- abs(mean(x[g == 1], na.rm = TRUE) - mean(x[g == 0], na.rm = TRUE))
  perms <- replicate(nperm, {
    gp <- sample(g)
    abs(mean(x[gp == 1], na.rm = TRUE) - mean(x[gp == 0], na.rm = TRUE))
  })
  mean(perms >= obs)
}

target_genes <- unique(c(selected_genes, late_genes, early_genes))

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
  member_map <- data.frame(
    geo_accession = sub("_.*$", "", members),
    member = members,
    stringsAsFactors = FALSE
  )
  meta_rif <- left_join(meta_rif, member_map, by = "geo_accession")

  out <- lapply(seq_len(nrow(meta_rif)), function(i) {
    mtx <- read_10x_from_nested_tar("path/to/local/external-storage/03_rif/GSE250130/GSE250130_RAW.tar", meta_rif$member[i])
    assign <- assign_cells_heca(mtx, centroid_mat)
    agg <- aggregate_target_counts(mtx, assign, target_genes)
    agg$sample <- meta_rif$sample[i]
    agg$group <- meta_rif$group_simple[i]
    agg$day <- meta_rif$day[i]
    cells <- assign %>% count(assigned_celltype, name = "cell_n") %>% mutate(sample = meta_rif$sample[i], group = meta_rif$group_simple[i])
    assign_sum <- assign %>% summarise(sample = meta_rif$sample[i], group = meta_rif$group_simple[i], mean_cor = mean(assignment_cor, na.rm = TRUE))
    list(agg = agg, cells = cells, sum = assign_sum)
  })
  list(
    agg = bind_rows(lapply(out, `[[`, "agg")),
    cells = bind_rows(lapply(out, `[[`, "cells")),
    summary = bind_rows(lapply(out, `[[`, "sum"))
  )
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
    agg <- aggregate_target_counts(mtx, assign, target_genes)
    agg$sample <- meta$geo_accession[i]
    agg$group <- meta$analysis_group[i]
    cells <- assign %>% count(assigned_celltype, name = "cell_n") %>% mutate(sample = meta$geo_accession[i], group = meta$analysis_group[i])
    assign_sum <- assign %>% summarise(sample = meta$geo_accession[i], group = meta$analysis_group[i], mean_cor = mean(assignment_cor, na.rm = TRUE))
    list(agg = agg, cells = cells, sum = assign_sum)
  })
  list(
    agg = bind_rows(lapply(out, `[[`, "agg")),
    cells = bind_rows(lapply(out, `[[`, "cells")),
    summary = bind_rows(lapply(out, `[[`, "sum"))
  )
}

process_endo213216 <- function() {
  score213 <- read.csv(file.path(root, "analysis/02_endometriosis/round4_deep/tables/GSE213216_signature_scores.csv"), check.names = FALSE)
  members <- utils::untar("path/to/local/external-storage/02_endometriosis/GSE213216/GSE213216_RAW.tar", list = TRUE)
  map213 <- data.frame(sample = sub("_.*$", "", members), member = members, stringsAsFactors = FALSE) %>%
    left_join(score213 %>% select(sample, category), by = "sample") %>%
    filter(!is.na(category))
  out <- lapply(seq_len(nrow(map213)), function(i) {
    mtx <- read_10x_from_nested_tar("path/to/local/external-storage/02_endometriosis/GSE213216/GSE213216_RAW.tar", map213$member[i], inner_prefix = sub("\\.tar\\.gz$", "", map213$member[i]))
    assign <- assign_cells_heca(mtx, centroid_mat)
    agg <- aggregate_target_counts(mtx, assign, target_genes)
    agg$sample <- map213$sample[i]
    agg$group <- map213$category[i]
    cells <- assign %>% count(assigned_celltype, name = "cell_n") %>% mutate(sample = map213$sample[i], group = map213$category[i])
    assign_sum <- assign %>% summarise(sample = map213$sample[i], group = map213$category[i], mean_cor = mean(assignment_cor, na.rm = TRUE))
    list(agg = agg, cells = cells, sum = assign_sum)
  })
  list(
    agg = bind_rows(lapply(out, `[[`, "agg")),
    cells = bind_rows(lapply(out, `[[`, "cells")),
    summary = bind_rows(lapply(out, `[[`, "sum"))
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
    agg <- aggregate_target_counts(mtx, assign, target_genes)
    agg$sample <- meta$geo_accession[i]
    agg$group <- meta$disease[i]
    agg$phase <- meta$phase[i]
    cells <- assign %>% count(assigned_celltype, name = "cell_n") %>% mutate(sample = meta$geo_accession[i], group = meta$disease[i], phase = meta$phase[i])
    assign_sum <- assign %>% summarise(sample = meta$geo_accession[i], group = meta$disease[i], phase = meta$phase[i], mean_cor = mean(assignment_cor, na.rm = TRUE))
    list(agg = agg, cells = cells, sum = assign_sum)
  })
  list(
    agg = bind_rows(lapply(out, `[[`, "agg")),
    cells = bind_rows(lapply(out, `[[`, "cells")),
    summary = bind_rows(lapply(out, `[[`, "sum"))
  )
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

rif <- load_or_build(file.path(tabdir, "rif_heca_raw_cache.rds"), process_rif)
endo179 <- load_or_build(file.path(tabdir, "endo179640_heca_raw_cache.rds"), process_endo179640)
endo213 <- load_or_build(file.path(tabdir, "endo213216_heca_raw_cache.rds"), process_endo213216)
endo214 <- load_or_build(file.path(tabdir, "endo214411_heca_raw_cache.rds"), process_endo214411)

write.csv(rif$cells, file.path(tabdir, "rif_heca_cellcounts.csv"), row.names = FALSE)
write.csv(endo179$cells, file.path(tabdir, "endo179640_heca_cellcounts.csv"), row.names = FALSE)
write.csv(endo213$cells, file.path(tabdir, "endo213216_heca_cellcounts.csv"), row.names = FALSE)
write.csv(endo214$cells, file.path(tabdir, "endo214411_heca_cellcounts.csv"), row.names = FALSE)

assign_summary <- bind_rows(
  rif$summary %>% mutate(dataset = "GSE250130"),
  endo179$summary %>% mutate(dataset = "GSE179640"),
  endo213$summary %>% mutate(dataset = "GSE213216"),
  endo214$summary %>% mutate(dataset = "GSE214411")
)
write.csv(assign_summary, file.path(tabdir, "heca_assignment_summary.csv"), row.names = FALSE)

rif_scores <- calc_scores_from_agg(rif$agg, c("sample", "group", "day"))
fertile_ref <- rif_scores %>% filter(!is.na(day))
celltype_models <- lapply(split(fertile_ref, fertile_ref$celltype), function(df) {
  if (nrow(df) >= 5 && length(unique(df$day)) >= 3) lm(day ~ receptivity_score, data = df) else NULL
})
global_model <- lm(day ~ receptivity_score, data = fertile_ref)
rif_scores <- bind_rows(lapply(split(rif_scores, rif_scores$celltype), function(df) {
  fit <- celltype_models[[unique(df$celltype)]]
  if (is.null(fit)) fit <- global_model
  proj <- project_to_timeline(df$receptivity_score, fit)
  bind_cols(df, proj)
}))
write.csv(rif_scores, file.path(tabdir, "rif_heca_celltype_timing_scores.csv"), row.names = FALSE)

project_external <- function(score_df) {
  bind_rows(lapply(split(score_df, score_df$celltype), function(df) {
    fit <- celltype_models[[unique(df$celltype)]]
    if (is.null(fit)) fit <- global_model
    proj <- project_to_timeline(df$receptivity_score, fit)
    bind_cols(df, proj)
  }))
}

endo179_scores <- calc_scores_from_agg(endo179$agg, c("sample", "group"))
endo179_scores <- project_external(endo179_scores)
write.csv(endo179_scores, file.path(tabdir, "endo179640_heca_celltype_timing_scores.csv"), row.names = FALSE)

endo213_scores <- calc_scores_from_agg(endo213$agg, c("sample", "group"))
endo213_scores <- project_external(endo213_scores)
write.csv(endo213_scores, file.path(tabdir, "endo213216_heca_celltype_timing_scores.csv"), row.names = FALSE)

endo214_scores <- calc_scores_from_agg(endo214$agg, c("sample", "group", "phase"))
endo214_scores <- project_external(endo214_scores)
write.csv(endo214_scores, file.path(tabdir, "endo214411_heca_celltype_timing_scores.csv"), row.names = FALSE)

# Figures
p_ref <- heca_stage %>%
  filter(stage_simple %in% c("Secretory_Early", "Secretory_EarlyMid", "Secretory_Mid", "Secretory_Late")) %>%
  mutate(stage_simple = factor(stage_simple, levels = c("Secretory_Early", "Secretory_EarlyMid", "Secretory_Mid", "Secretory_Late"))) %>%
  ggplot(aes(stage_simple, broad_celltype, fill = n_cells)) +
  geom_tile() +
  scale_fill_gradient(low = "#F4E8C1", high = "#A73A24") +
  theme_bw(base_size = 11) +
  labs(title = "HECA secretory-stage coverage across broad endometrial cell types", x = "", y = "")
ggsave(file.path(figdir, "Figure_ATLAS_HECA_1_reference_stage_coverage.png"), p_ref, width = 8.2, height = 4.8, dpi = 300)

p_rif <- rif_scores %>%
  filter(group %in% c("Fertile_LH7", "RIF"), celltype != "Unknown") %>%
  ggplot(aes(celltype, pred_day, color = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75), alpha = 0.7, size = 1.15) +
  geom_hline(yintercept = 7, linetype = 2, color = "grey45") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "HECA-anchored cell-type-resolved timing projection in RIF", x = "", y = "Predicted receptive day")
ggsave(file.path(figdir, "Figure_ATLAS_HECA_2_rif_celltype_timing.png"), p_rif, width = 10.6, height = 5.6, dpi = 300)

endo_combined <- bind_rows(
  endo179_scores %>% mutate(dataset = "GSE179640"),
  endo213_scores %>% mutate(dataset = "GSE213216"),
  endo214_scores %>% mutate(dataset = "GSE214411")
)
p_endo <- endo_combined %>%
  filter(celltype %in% c("Luminal_Epi", "Glandular_Epi", "Stroma", "Decidual_Stroma", "Lymphoid", "Myeloid", "Endothelial")) %>%
  mutate(group2 = case_when(
    group %in% c("Control_Eutopic", "Control") ~ "Control",
    group %in% c("Endo_Eutopic", "Eutopic", "endometriosis") ~ "Eutopic_endo",
    group %in% c("Ectopic", "Lesion") ~ "Lesion",
    group %in% c("Ectopic_Ovary", "Ovary") ~ "Ovary",
    TRUE ~ as.character(group)
  )) %>%
  filter(group2 %in% c("Control", "Eutopic_endo", "Lesion", "Ovary")) %>%
  ggplot(aes(celltype, pred_day, color = group2)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.78)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.78), alpha = 0.65, size = 0.9) +
  facet_wrap(~ dataset, ncol = 1) +
  geom_hline(yintercept = 7, linetype = 2, color = "grey45") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(title = "HECA-anchored cell-type timing across endometriosis cohorts", x = "", y = "Predicted receptive day")
ggsave(file.path(figdir, "Figure_ATLAS_HECA_3_endo_celltype_timing.png"), p_endo, width = 10.8, height = 11.0, dpi = 300)

comp_df <- bind_rows(
  rif$cells %>% mutate(dataset = "GSE250130"),
  endo179$cells %>% mutate(dataset = "GSE179640"),
  endo213$cells %>% mutate(dataset = "GSE213216"),
  endo214$cells %>% mutate(dataset = "GSE214411")
) %>%
  group_by(dataset, sample) %>%
  mutate(prop = cell_n / sum(cell_n)) %>%
  ungroup()

p_comp <- comp_df %>%
  filter(assigned_celltype != "Unknown") %>%
  ggplot(aes(sample, prop, fill = assigned_celltype)) +
  geom_col() +
  facet_wrap(~ dataset, scales = "free_x", ncol = 1) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "HECA-based broad cell-type composition across the four key single-cell cohorts", x = "Samples", y = "Cell fraction", fill = "Cell type")
ggsave(file.path(figdir, "Figure_ATLAS_HECA_4_celltype_composition.png"), p_comp, width = 10.8, height = 9.8, dpi = 300)

heat_df <- bind_rows(
  rif_scores %>% mutate(dataset = "RIF") %>% transmute(label = paste(dataset, group, sep = " | "), celltype, pred_day),
  endo179_scores %>% mutate(dataset = "Endo179640") %>% transmute(label = paste(dataset, group, sep = " | "), celltype, pred_day),
  endo213_scores %>% mutate(dataset = "Endo213216") %>% transmute(label = paste(dataset, group, sep = " | "), celltype, pred_day),
  endo214_scores %>% mutate(dataset = "Endo214411") %>% transmute(label = paste(dataset, group, sep = " | "), celltype, pred_day)
) %>%
  group_by(label, celltype) %>%
  summarise(mean_pred_day = mean(pred_day, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = label, values_from = mean_pred_day)
heat_mat <- as.data.frame(heat_df)
rownames(heat_mat) <- heat_mat$celltype
heat_mat$celltype <- NULL
png(file.path(figdir, "Figure_ATLAS_HECA_5_cross_cohort_heatmap.png"), width = 2600, height = 1700, res = 280)
pheatmap(as.matrix(heat_mat), cluster_rows = FALSE, cluster_cols = FALSE, main = "Mean predicted receptive day by HECA-mapped cell type")
dev.off()

# Statistics
rif_stats <- rif_scores %>%
  filter(group %in% c("Fertile_LH7", "RIF"), celltype != "Unknown") %>%
  group_by(celltype) %>%
  summarise(
    estimate = bootstrap_diff(pred_day[group == "RIF"], pred_day[group == "Fertile_LH7"])[["obs"]],
    lower = bootstrap_diff(pred_day[group == "RIF"], pred_day[group == "Fertile_LH7"])[["lower"]],
    upper = bootstrap_diff(pred_day[group == "RIF"], pred_day[group == "Fertile_LH7"])[["upper"]],
    p_perm = perm_p(pred_day, ifelse(group == "RIF", 1, 0)),
    .groups = "drop"
  )
write.csv(rif_stats, file.path(tabdir, "rif_heca_celltype_stats.csv"), row.names = FALSE)

endo_stats <- bind_rows(
  endo179_scores %>% filter(group %in% c("Control_Eutopic", "Endo_Eutopic")) %>%
    group_by(celltype) %>%
    summarise(
      dataset = "GSE179640",
      comparison = "Endo_Eutopic vs Control_Eutopic",
      estimate = bootstrap_diff(pred_day[group == "Endo_Eutopic"], pred_day[group == "Control_Eutopic"])[["obs"]],
      lower = bootstrap_diff(pred_day[group == "Endo_Eutopic"], pred_day[group == "Control_Eutopic"])[["lower"]],
      upper = bootstrap_diff(pred_day[group == "Endo_Eutopic"], pred_day[group == "Control_Eutopic"])[["upper"]],
      p_perm = perm_p(pred_day, ifelse(group == "Endo_Eutopic", 1, 0)),
      .groups = "drop"
    ),
  endo214_scores %>% filter(group %in% c("Control", "endometriosis")) %>%
    group_by(celltype) %>%
    summarise(
      dataset = "GSE214411",
      comparison = "Endometriosis vs Control",
      estimate = bootstrap_diff(pred_day[group == "endometriosis"], pred_day[group == "Control"])[["obs"]],
      lower = bootstrap_diff(pred_day[group == "endometriosis"], pred_day[group == "Control"])[["lower"]],
      upper = bootstrap_diff(pred_day[group == "endometriosis"], pred_day[group == "Control"])[["upper"]],
      p_perm = perm_p(pred_day, ifelse(group == "endometriosis", 1, 0)),
      .groups = "drop"
    )
)
write.csv(endo_stats, file.path(tabdir, "endo_heca_celltype_stats.csv"), row.names = FALSE)

summary_tbl <- bind_rows(
  rif_stats %>% mutate(dataset = "GSE250130", comparison = "RIF vs Fertile_LH7"),
  endo_stats
) %>%
  select(dataset, comparison, celltype, estimate, lower, upper, p_perm)
write.csv(summary_tbl, file.path(tabdir, "heca_key_celltype_shift_summary.csv"), row.names = FALSE)

coverage_tbl <- bind_rows(
  data.frame(dataset = "GSE250130", metric = "Mean assignment correlation", value = mean(rif$summary$mean_cor, na.rm = TRUE)),
  data.frame(dataset = "GSE179640", metric = "Mean assignment correlation", value = mean(endo179$summary$mean_cor, na.rm = TRUE)),
  data.frame(dataset = "GSE213216", metric = "Mean assignment correlation", value = mean(endo213$summary$mean_cor, na.rm = TRUE)),
  data.frame(dataset = "GSE214411", metric = "Mean assignment correlation", value = mean(endo214$summary$mean_cor, na.rm = TRUE))
)
write.csv(coverage_tbl, file.path(tabdir, "heca_assignment_quality_summary.csv"), row.names = FALSE)