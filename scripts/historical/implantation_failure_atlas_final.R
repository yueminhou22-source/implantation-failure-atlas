#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(Matrix)
  library(jsonlite)
  library(pheatmap)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(hugene10sttranscriptcluster.db)
})

outdir <- "analysis/05_implantation_failure_atlas/final"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

compute_signature_score <- function(expr_mat, up_genes, down_genes = character()) {
  up <- intersect(up_genes, rownames(expr_mat))
  down <- intersect(down_genes, rownames(expr_mat))
  if (length(up) == 0 && length(down) == 0) return(rep(NA_real_, ncol(expr_mat)))
  z <- t(scale(t(expr_mat)))
  z[!is.finite(z)] <- NA_real_
  up_score <- if (length(up) > 0) colMeans(z[up, , drop = FALSE], na.rm = TRUE) else rep(0, ncol(expr_mat))
  down_score <- if (length(down) > 0) colMeans(z[down, , drop = FALSE], na.rm = TRUE) else rep(0, ncol(expr_mat))
  up_score - down_score
}

compute_signature_score_raw <- function(expr_mat, up_genes, down_genes = character()) {
  up <- intersect(up_genes, rownames(expr_mat))
  down <- intersect(down_genes, rownames(expr_mat))
  if (length(up) == 0 && length(down) == 0) return(rep(NA_real_, ncol(expr_mat)))
  up_score <- if (length(up) > 0) colMeans(expr_mat[up, , drop = FALSE], na.rm = TRUE) else rep(0, ncol(expr_mat))
  down_score <- if (length(down) > 0) colMeans(expr_mat[down, , drop = FALSE], na.rm = TRUE) else rep(0, ncol(expr_mat))
  up_score - down_score
}

read_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  table_start <- which(lines == "!series_matrix_table_begin")
  table_end <- which(lines == "!series_matrix_table_end")
  if (length(table_start) == 0 || length(table_end) == 0) stop("Could not find series matrix table in ", path)
  table_txt <- paste(lines[(table_start + 1):(table_end - 1)], collapse = "\n")
  mat <- read.delim(text = table_txt, check.names = FALSE, stringsAsFactors = FALSE)
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
  list(expr = mat, meta = meta)
}

read_series_meta <- function(path) {
  read_series_matrix(path)$meta
}

collapse_entrez_expr <- function(expr_df) {
  names(expr_df)[1] <- "entrez"
  expr_df$entrez <- as.character(expr_df$entrez)
  expr_df$symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = expr_df$entrez,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  expr_df <- expr_df %>% filter(!is.na(symbol), symbol != "")
  expr_df <- expr_df %>%
    dplyr::select(-entrez) %>%
    group_by(symbol) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  expr_mat <- as.data.frame(expr_df)
  rownames(expr_mat) <- expr_mat$symbol
  expr_mat$symbol <- NULL
  expr_mat
}

extract_gse111974_probe_map <- function(raw_tar) {
  exdir <- tempfile("gse111974_")
  dir.create(exdir)
  utils::untar(raw_tar, exdir = exdir, files = "GSM3045867_SG12324209_253949426931_S001_GE1_1100_Jul11_1_1.txt.gz")
  raw_path <- file.path(exdir, "GSM3045867_SG12324209_253949426931_S001_GE1_1100_Jul11_1_1.txt.gz")
  raw_lines <- readLines(gzfile(raw_path), warn = FALSE)
  feat_idx <- grep("^FEATURES", raw_lines)[1]
  data_idx <- grep("^DATA", raw_lines)
  data_idx <- data_idx[data_idx > feat_idx]
  header <- strsplit(raw_lines[feat_idx], "\t")[[1]]
  feature_lines <- raw_lines[data_idx]
  feat_txt <- paste(c(paste(header, collapse = "\t"), feature_lines), collapse = "\n")
  feat_df <- read.delim(text = feat_txt, check.names = FALSE, stringsAsFactors = FALSE)
  feat_df %>%
    filter(ControlType == 0) %>%
    transmute(probe_id = ProbeName, refseq = SystematicName) %>%
    distinct()
}

collapse_expr_by_gene <- function(expr_df, probe_map, probe_col = "ID_REF") {
  names(expr_df)[1] <- probe_col
  merged <- inner_join(probe_map, expr_df, by = setNames(probe_col, "probe_id"))
  merged$symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = merged$refseq,
    keytype = "REFSEQ",
    column = "SYMBOL",
    multiVals = "first"
  )
  merged <- merged %>% filter(!is.na(symbol), symbol != "")
  expr_only <- merged %>% select(-refseq, -probe_id)
  expr_collapsed <- expr_only %>%
    group_by(symbol) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  expr_mat <- as.data.frame(expr_collapsed)
  rownames(expr_mat) <- expr_mat$symbol
  expr_mat$symbol <- NULL
  expr_mat
}

extract_simple_sample_pseudobulk <- function(raw_tar, member_prefix) {
  exdir <- tempfile("gse214411_")
  dir.create(exdir)
  members <- c(
    paste0(member_prefix, "_features.tsv.gz"),
    paste0(member_prefix, "_barcodes.tsv.gz"),
    paste0(member_prefix, "_matrix.mtx.gz")
  )
  utils::untar(raw_tar, exdir = exdir, files = members)
  mtx <- readMM(file.path(exdir, paste0(member_prefix, "_matrix.mtx.gz")))
  features <- read.delim(gzfile(file.path(exdir, paste0(member_prefix, "_features.tsv.gz"))), header = FALSE, stringsAsFactors = FALSE)
  rownames(mtx) <- make.unique(features[[2]])
  counts <- Matrix::rowSums(mtx)
  data.frame(gene = names(counts), value = as.numeric(counts), stringsAsFactors = FALSE)
}

parse_geomx_counts <- function(raw_tar, pkc_gz) {
  pkc <- jsonlite::fromJSON(gzfile(pkc_gz))
  probe_map <- bind_rows(lapply(seq_len(nrow(pkc$Targets)), function(i) {
    pr <- pkc$Targets$Probes[[i]]
    if (is.null(pr) || nrow(pr) == 0) return(NULL)
    data.frame(
      rts_id = pr$RTS_ID,
      symbol = rep(pkc$Targets$DisplayName[i], nrow(pr)),
      stringsAsFactors = FALSE
    )
  })) %>% filter(!is.na(symbol), symbol != "")
  members <- utils::untar(raw_tar, list = TRUE)
  dcc_members <- members[grepl("\\.dcc\\.gz$", members)]
  score_list <- lapply(dcc_members, function(member) {
    exdir <- tempfile("geomx_")
    dir.create(exdir)
    utils::untar(raw_tar, exdir = exdir, files = member)
    dcc_path <- file.path(exdir, member)
    lines <- readLines(gzfile(dcc_path), warn = FALSE)
    start <- grep("^<Code_Summary>", lines)[1]
    end <- grep("^</Code_Summary>", lines)[1]
    dat <- read.csv(text = paste(lines[(start + 1):(end - 1)], collapse = "\n"), header = FALSE, stringsAsFactors = FALSE)
    colnames(dat) <- c("rts_id", "count")
    dat <- inner_join(probe_map, dat, by = "rts_id") %>%
      group_by(symbol) %>%
      summarise(count = sum(as.numeric(count), na.rm = TRUE), .groups = "drop")
    gsm <- sub("_.*$", "", basename(member))
    data.frame(gene = dat$symbol, value = dat$count, sample = gsm, stringsAsFactors = FALSE)
  })
  bind_rows(score_list)
}

compute_effect_size <- function(score, group) {
  idx1 <- which(group == 1)
  idx0 <- which(group == 0)
  m1 <- mean(score[idx1], na.rm = TRUE)
  m0 <- mean(score[idx0], na.rm = TRUE)
  s1 <- sd(score[idx1], na.rm = TRUE)
  s0 <- sd(score[idx0], na.rm = TRUE)
  n1 <- length(idx1)
  n0 <- length(idx0)
  sp <- sqrt(((n1 - 1) * s1^2 + (n0 - 1) * s0^2) / pmax(n1 + n0 - 2, 1))
  d <- (m1 - m0) / sp
  se <- sqrt((n1 + n0) / (n1 * n0) + (d^2 / (2 * pmax(n1 + n0 - 2, 1))))
  data.frame(effect = d, se = se, lower = d - 1.96 * se, upper = d + 1.96 * se, n_case = n1, n_ctrl = n0)
}

random_effects_meta <- function(df) {
  w <- 1 / (df$se^2)
  fixed <- sum(w * df$effect) / sum(w)
  q <- sum(w * (df$effect - fixed)^2)
  cval <- sum(w) - sum(w^2) / sum(w)
  tau2 <- max((q - (nrow(df) - 1)) / cval, 0)
  w_re <- 1 / (df$se^2 + tau2)
  pooled <- sum(w_re * df$effect) / sum(w_re)
  se_pooled <- sqrt(1 / sum(w_re))
  data.frame(effect = pooled, se = se_pooled, lower = pooled - 1.96 * se_pooled, upper = pooled + 1.96 * se_pooled, tau2 = tau2)
}

project_to_timeline <- function(score, calib_model) {
  pred <- predict(calib_model, newdata = data.frame(receptivity_score = score))
  state <- cut(pred, breaks = c(-Inf, 6.5, 8.5, Inf), labels = c("Delayed", "In-phase", "Advanced"))
  data.frame(pred_day = pred, woi_distance = abs(pred - 7), timing_state = state)
}

decidualization_genes <- c("IGFBP1", "PRL", "LEFTY2", "FOXO1", "WNT4", "HAND2", "IL15", "SPP1", "PAEP", "GPX3")
hormone_response_genes <- c("PGR", "ESR1", "GREB1", "IHH", "HOXA10", "HOXA11", "KLF9", "NR2F2", "MUC1", "CXCL14")

receptivity_sig <- read.csv("analysis/03_rif/round5_maximal/tables/receptivity_signature_genes_round5.csv", check.names = FALSE)
late_genes <- na.omit(receptivity_sig$late_genes)
early_genes <- na.omit(receptivity_sig$early_genes)
rif_sig <- read.csv("analysis/03_rif/round5_maximal/tables/rif_directional_signature_genes_round5.csv", check.names = FALSE)
rif_up_genes <- na.omit(rif_sig$rif_up)
rif_down_genes <- na.omit(rif_sig$rif_down)
adeno_cons <- read.csv("analysis/04_adenomyosis/round4_deep/tables/consensus_adenomyosis_signature.csv", check.names = FALSE)
adeno_up <- na.omit(adeno_cons$consensus_up)
adeno_down <- na.omit(adeno_cons$consensus_down)
lesion_deg <- read.csv("analysis/02_endometriosis/round4_deep/tables/Figure_ENDO_2A_Ectopic_vs_Eutopic_deg.csv", check.names = FALSE)
lesion_up <- lesion_deg %>% filter(adj.P.Val < 0.05, logFC > 0) %>% arrange(adj.P.Val) %>% slice_head(n = 30) %>% pull(gene)
lesion_down <- lesion_deg %>% filter(adj.P.Val < 0.05, logFC < 0) %>% arrange(adj.P.Val) %>% slice_head(n = 30) %>% pull(gene)

score_rif <- read.csv("analysis/03_rif/round5_maximal/tables/receptivity_score_by_sample_round5.csv", check.names = FALSE)
calib_model <- lm(day ~ receptivity_score, data = score_rif %>% filter(!is.na(day)))
rif_proj <- read.csv("analysis/03_rif/round5_maximal/tables/rif_projection_to_fertile_timeline_round5.csv", check.names = FALSE)

all_scores <- list()

# RIF discovery and external
rif_discovery <- rif_proj %>%
  transmute(
    dataset = "GSE250130",
    disease = "RIF",
    cohort_type = "Single-cell pseudobulk",
    sample = sample,
    group = group_simple,
    tissue_context = ifelse(group_simple == "RIF", "RIF", "Fertile_reference"),
    receptivity_score = receptivity_score,
    pred_day = pred_day,
    woi_distance = abs(pred_day - 7),
    timing_state = as.character(timing_state),
    rif_immune_score = NA_real_,
    decidualization_score = NA_real_,
    hormone_score = NA_real_,
    lesion_score = NA_real_,
    adeno_consensus_score = NA_real_
  )
all_scores[["rif_discovery"]] <- rif_discovery

score111 <- read.csv("analysis/03_rif/round5_maximal/tables/GSE111974_signature_scores_round5.csv", check.names = FALSE)
proj111 <- project_to_timeline(score111$receptivity_score, calib_model)
all_scores[["rif_gse111974"]] <- bind_cols(
  data.frame(
    dataset = "GSE111974", disease = "RIF", cohort_type = "Bulk endometrium",
    sample = score111$sample, group = score111$group, tissue_context = score111$group,
    receptivity_score = score111$receptivity_score,
    timing_state = as.character(proj111$timing_state),
    rif_immune_score = score111$rif_immune_score,
    decidualization_score = NA_real_, hormone_score = NA_real_,
    lesion_score = NA_real_, adeno_consensus_score = NA_real_
  ),
  proj111[, c("pred_day", "woi_distance")]
)

score58144 <- read.csv("analysis/03_rif/round5_maximal/tables/GSE58144_signature_scores_round5.csv", check.names = FALSE)
proj58144 <- project_to_timeline(score58144$receptivity_score, calib_model)
all_scores[["rif_gse58144"]] <- bind_cols(
  data.frame(
    dataset = "GSE58144", disease = "RIF", cohort_type = "Bulk endometrium",
    sample = score58144$sample, group = score58144$group, tissue_context = score58144$group,
    receptivity_score = score58144$receptivity_score,
    timing_state = as.character(proj58144$timing_state),
    rif_immune_score = score58144$rif_immune_score,
    decidualization_score = NA_real_, hormone_score = NA_real_,
    lesion_score = NA_real_, adeno_consensus_score = NA_real_
  ),
  proj58144[, c("pred_day", "woi_distance")]
)

score_rif_spatial <- read.csv("analysis/03_rif/round5_maximal/tables/GSE287278_spatial_signature_scores_round5.csv", check.names = FALSE)
proj_rif_spatial <- project_to_timeline(score_rif_spatial$receptivity_score, calib_model)
all_scores[["rif_spatial"]] <- bind_cols(
  data.frame(
    dataset = "GSE287278", disease = "RIF", cohort_type = "Spatial transcriptomics",
    sample = paste(score_rif_spatial$sample, score_rif_spatial$barcode, sep = "::"),
    group = score_rif_spatial$group, tissue_context = score_rif_spatial$group,
    receptivity_score = score_rif_spatial$receptivity_score,
    timing_state = as.character(proj_rif_spatial$timing_state),
    rif_immune_score = score_rif_spatial$rif_immune_score,
    decidualization_score = NA_real_, hormone_score = NA_real_,
    lesion_score = NA_real_, adeno_consensus_score = NA_real_
  ),
  proj_rif_spatial[, c("pred_day", "woi_distance")]
)

# Endometriosis discovery pseudobulk
endo_counts <- read.csv("analysis/02_endometriosis/tables/GSE179640_firstpass/pseudobulk_counts.csv", check.names = FALSE)
colnames(endo_counts)[1] <- "gene"
endo_counts <- endo_counts %>% group_by(gene) %>% summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
endo_expr <- as.data.frame(endo_counts)
rownames(endo_expr) <- endo_expr$gene
endo_expr$gene <- NULL
endo_expr <- log2(as.matrix(endo_expr) + 1)
endo_meta <- read.csv("analysis/02_endometriosis/metadata/GSE179640_metadata_round2.csv", check.names = FALSE)
endo_score <- data.frame(
  sample = colnames(endo_expr),
  group = endo_meta$analysis_group[match(colnames(endo_expr), endo_meta$sample)],
  receptivity_score = compute_signature_score(endo_expr, late_genes, early_genes),
  rif_immune_score = compute_signature_score(endo_expr, rif_up_genes, rif_down_genes),
  decidualization_score = compute_signature_score(endo_expr, decidualization_genes, character()),
  hormone_score = compute_signature_score(endo_expr, hormone_response_genes, character()),
  lesion_score = compute_signature_score(endo_expr, lesion_up, lesion_down),
  adeno_consensus_score = compute_signature_score(endo_expr, adeno_up, adeno_down),
  stringsAsFactors = FALSE
)
proj_endo <- project_to_timeline(endo_score$receptivity_score, calib_model)
all_scores[["endo_gse179640"]] <- bind_cols(
  data.frame(
    dataset = "GSE179640", disease = "Endometriosis", cohort_type = "Single-cell pseudobulk",
    sample = endo_score$sample, group = endo_score$group, tissue_context = endo_score$group,
    receptivity_score = endo_score$receptivity_score,
    timing_state = as.character(proj_endo$timing_state),
    rif_immune_score = endo_score$rif_immune_score,
    decidualization_score = endo_score$decidualization_score,
    hormone_score = endo_score$hormone_score,
    lesion_score = endo_score$lesion_score,
    adeno_consensus_score = endo_score$adeno_consensus_score
  ),
  proj_endo[, c("pred_day", "woi_distance")]
)

# Endometriosis independent sc pseudobulk
g213 <- read.csv("analysis/02_endometriosis/round4_deep/tables/GSE213216_pseudobulk_counts.csv", check.names = FALSE)
colnames(g213)[1] <- "gene"
g213_expr <- as.data.frame(g213)
rownames(g213_expr) <- g213_expr$gene
g213_expr$gene <- NULL
g213_expr <- log2(as.matrix(g213_expr) + 1)
g213_meta <- read.csv("analysis/02_endometriosis/round4_deep/tables/GSE213216_signature_scores.csv", check.names = FALSE)
g213_score <- data.frame(
  sample = colnames(g213_expr),
  group = g213_meta$category[match(colnames(g213_expr), g213_meta$sample)],
  receptivity_score = compute_signature_score(g213_expr, late_genes, early_genes),
  rif_immune_score = compute_signature_score(g213_expr, rif_up_genes, rif_down_genes),
  decidualization_score = compute_signature_score(g213_expr, decidualization_genes, character()),
  hormone_score = compute_signature_score(g213_expr, hormone_response_genes, character()),
  lesion_score = compute_signature_score(g213_expr, lesion_up, lesion_down),
  adeno_consensus_score = compute_signature_score(g213_expr, adeno_up, adeno_down),
  stringsAsFactors = FALSE
)
proj_g213 <- project_to_timeline(g213_score$receptivity_score, calib_model)
all_scores[["endo_gse213216"]] <- bind_cols(
  data.frame(
    dataset = "GSE213216", disease = "Endometriosis", cohort_type = "Single-cell pseudobulk",
    sample = g213_score$sample, group = g213_score$group, tissue_context = g213_score$group,
    receptivity_score = g213_score$receptivity_score,
    timing_state = as.character(proj_g213$timing_state),
    rif_immune_score = g213_score$rif_immune_score,
    decidualization_score = g213_score$decidualization_score,
    hormone_score = g213_score$hormone_score,
    lesion_score = g213_score$lesion_score,
    adeno_consensus_score = g213_score$adeno_consensus_score
  ),
  proj_g213[, c("pred_day", "woi_distance")]
)

# Endometriosis eutopic validation GSE214411
meta214a <- read_series_meta("/Volumes/Extreme SSD/02_endometriosis/GSE214411/GSE214411-GPL24676_series_matrix.txt.gz")
meta214a$disease <- sub("^disease: ", "", meta214a$c2)
meta214a$phase <- sub("^cycle phase: ", "", meta214a$c3)
meta214b <- read_series_meta("/Volumes/Extreme SSD/02_endometriosis/GSE214411/GSE214411-GPL11154_series_matrix.txt.gz")
meta214b$disease <- "Control"
meta214b$phase <- NA_character_
meta214 <- bind_rows(meta214a, meta214b)
tar_members214 <- utils::untar("/Volumes/Extreme SSD/02_endometriosis/GSE263897/GSE214411_RAW.tar", list = TRUE)
member_prefixes214 <- unique(sub("_(features|barcodes|matrix)\\.tsv\\.gz$|_matrix\\.mtx\\.gz$", "", tar_members214))
member_df214 <- data.frame(member_prefix = member_prefixes214, stringsAsFactors = FALSE)
member_df214$geo_accession <- sub("_.*$", "", member_df214$member_prefix)
meta214 <- left_join(meta214, member_df214, by = "geo_accession")
pb214_list <- lapply(meta214$member_prefix, function(stem) extract_simple_sample_pseudobulk("/Volumes/Extreme SSD/02_endometriosis/GSE263897/GSE214411_RAW.tar", stem))
for (i in seq_along(pb214_list)) names(pb214_list[[i]])[2] <- meta214$geo_accession[i]
pb214 <- Reduce(function(x, y) full_join(x, y, by = "gene"), pb214_list)
pb214[is.na(pb214)] <- 0
pb214_mat <- as.data.frame(pb214)
rownames(pb214_mat) <- pb214_mat$gene
pb214_mat$gene <- NULL
pb214_mat <- log2(as.matrix(pb214_mat) + 1)
score214 <- data.frame(
  sample = colnames(pb214_mat),
  group = meta214$disease[match(colnames(pb214_mat), meta214$geo_accession)],
  tissue_context = meta214$phase[match(colnames(pb214_mat), meta214$geo_accession)],
  receptivity_score = compute_signature_score(pb214_mat, late_genes, early_genes),
  rif_immune_score = compute_signature_score(pb214_mat, rif_up_genes, rif_down_genes),
  decidualization_score = compute_signature_score(pb214_mat, decidualization_genes, character()),
  hormone_score = compute_signature_score(pb214_mat, hormone_response_genes, character()),
  lesion_score = compute_signature_score(pb214_mat, lesion_up, lesion_down),
  adeno_consensus_score = compute_signature_score(pb214_mat, adeno_up, adeno_down),
  stringsAsFactors = FALSE
)
proj214 <- project_to_timeline(score214$receptivity_score, calib_model)
all_scores[["endo_gse214411"]] <- bind_cols(
  data.frame(
    dataset = "GSE214411", disease = "Endometriosis", cohort_type = "Single-cell pseudobulk",
    sample = score214$sample, group = score214$group, tissue_context = ifelse(is.na(score214$tissue_context), score214$group, score214$tissue_context),
    receptivity_score = score214$receptivity_score,
    timing_state = as.character(proj214$timing_state),
    rif_immune_score = score214$rif_immune_score,
    decidualization_score = score214$decidualization_score,
    hormone_score = score214$hormone_score,
    lesion_score = score214$lesion_score,
    adeno_consensus_score = score214$adeno_consensus_score
  ),
  proj214[, c("pred_day", "woi_distance")]
)

# Endometriosis bulk GSE135485
meta135 <- read_series_meta("/Volumes/Extreme SSD/02_endometriosis/GSE135485/GSE135485_series_matrix.txt.gz")
meta135$status <- sub("^subject status: ", "", meta135$c1)
expr135 <- read.csv(gzfile("/Volumes/Extreme SSD/02_endometriosis/GSE135485/GSE135485_Endometriosis_raw_counts.csv.gz"), check.names = FALSE)
names(expr135)[1] <- "gene"
expr135 <- expr135 %>% group_by(gene) %>% summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
expr135_mat <- as.data.frame(expr135)
rownames(expr135_mat) <- expr135_mat$gene
expr135_mat$gene <- NULL
expr135_mat <- log2(as.matrix(expr135_mat) + 1)
score135 <- data.frame(
  sample = colnames(expr135_mat),
  group = meta135$status[match(colnames(expr135_mat), meta135$title)],
  receptivity_score = compute_signature_score(expr135_mat, late_genes, early_genes),
  rif_immune_score = compute_signature_score(expr135_mat, rif_up_genes, rif_down_genes),
  decidualization_score = compute_signature_score(expr135_mat, decidualization_genes, character()),
  hormone_score = compute_signature_score(expr135_mat, hormone_response_genes, character()),
  lesion_score = compute_signature_score(expr135_mat, lesion_up, lesion_down),
  adeno_consensus_score = compute_signature_score(expr135_mat, adeno_up, adeno_down),
  stringsAsFactors = FALSE
)
proj135 <- project_to_timeline(score135$receptivity_score, calib_model)
all_scores[["endo_gse135485"]] <- bind_cols(
  data.frame(
    dataset = "GSE135485", disease = "Endometriosis", cohort_type = "Bulk endometrium",
    sample = score135$sample, group = score135$group, tissue_context = score135$group,
    receptivity_score = score135$receptivity_score,
    timing_state = as.character(proj135$timing_state),
    rif_immune_score = score135$rif_immune_score,
    decidualization_score = score135$decidualization_score,
    hormone_score = score135$hormone_score,
    lesion_score = score135$lesion_score,
    adeno_consensus_score = score135$adeno_consensus_score
  ),
  proj135[, c("pred_day", "woi_distance")]
)

# Endometriosis spatial GSE263897
meta263 <- read_series_meta("/Volumes/Extreme SSD/02_endometriosis/GSE263897/GSE263897_series_matrix.txt.gz")
meta263$tissue <- sub("^tissue: ", "", meta263$c1)
meta263$sample_id <- sub("^sampleID: ", "", meta263$c2)
meta263$cell_type <- sub("^cell type: ", "", meta263$c3)
geomx_long <- parse_geomx_counts("/Volumes/Extreme SSD/02_endometriosis/GSE263897/GSE263897_RAW.tar",
                                 "/Volumes/Extreme SSD/02_endometriosis/GSE263897/GSE263897_Hs_R_NGS_WTA_v1.0.pkc.gz")
geomx_mat <- geomx_long %>% pivot_wider(names_from = sample, values_from = value)
geomx_mat <- as.data.frame(geomx_mat)
rownames(geomx_mat) <- geomx_mat$gene
geomx_mat$gene <- NULL
geomx_mat[is.na(geomx_mat)] <- 0
geomx_mat <- log2(as.matrix(geomx_mat) + 1)
score263 <- data.frame(
  sample = colnames(geomx_mat),
  group = meta263$tissue[match(colnames(geomx_mat), meta263$geo_accession)],
  tissue_context = meta263$cell_type[match(colnames(geomx_mat), meta263$geo_accession)],
  receptivity_score = compute_signature_score(geomx_mat, late_genes, early_genes),
  rif_immune_score = compute_signature_score(geomx_mat, rif_up_genes, rif_down_genes),
  decidualization_score = compute_signature_score(geomx_mat, decidualization_genes, character()),
  hormone_score = compute_signature_score(geomx_mat, hormone_response_genes, character()),
  lesion_score = compute_signature_score(geomx_mat, lesion_up, lesion_down),
  adeno_consensus_score = compute_signature_score(geomx_mat, adeno_up, adeno_down),
  stringsAsFactors = FALSE
)
proj263 <- project_to_timeline(score263$receptivity_score, calib_model)
all_scores[["endo_gse263897"]] <- bind_cols(
  data.frame(
    dataset = "GSE263897", disease = "Endometriosis", cohort_type = "Spatial transcriptomics",
    sample = score263$sample, group = score263$group, tissue_context = score263$tissue_context,
    receptivity_score = score263$receptivity_score,
    timing_state = as.character(proj263$timing_state),
    rif_immune_score = score263$rif_immune_score,
    decidualization_score = score263$decidualization_score,
    hormone_score = score263$hormone_score,
    lesion_score = score263$lesion_score,
    adeno_consensus_score = score263$adeno_consensus_score
  ),
  proj263[, c("pred_day", "woi_distance")]
)

# Adenomyosis organoid discovery
expr244_raw <- as.data.frame(read_excel("/Volumes/Extreme SSD/04_adenomyosis/GSE244236/GSE244236_Normalized_counts.xlsx"))
gene_ids244 <- as.character(expr244_raw[[1]])
expr244 <- expr244_raw[, -1, drop = FALSE]
expr244[] <- lapply(expr244, as.numeric)
rownames(expr244) <- gene_ids244
id_map244 <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(expr244), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
expr244$symbol <- unname(id_map244)
expr244 <- expr244 %>% filter(!is.na(symbol), symbol != "")
expr2442 <- expr244 %>% dplyr::select(symbol, everything()) %>% group_by(symbol) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
expr244_mat <- as.data.frame(expr2442)
rownames(expr244_mat) <- expr244_mat$symbol
expr244_mat$symbol <- NULL
expr244_mat <- as.matrix(expr244_mat)
mode(expr244_mat) <- "numeric"
meta244 <- read.csv("analysis/04_adenomyosis/metadata/GSE244236_metadata.csv", check.names = FALSE)
meta244$sample_id <- sub(".*\\[([^]]+)\\].*", "\\1", meta244$title)
meta244$group <- ifelse(grepl("control", meta244$title, ignore.case = TRUE), "Control", "Adenomyosis")
meta244$phase_short <- ifelse(meta244$differentiation_phase == "mid-secretory phase", "SEC", "GEST")
common244 <- intersect(colnames(expr244_mat), meta244$sample_id)
expr244_mat <- expr244_mat[, common244, drop = FALSE]
meta244 <- meta244[match(common244, meta244$sample_id), ]
score244 <- data.frame(
  sample = colnames(expr244_mat),
  group = meta244$group,
  tissue_context = meta244$phase_short,
  receptivity_score = compute_signature_score(expr244_mat, late_genes, early_genes),
  rif_immune_score = compute_signature_score(expr244_mat, rif_up_genes, rif_down_genes),
  decidualization_score = compute_signature_score(expr244_mat, decidualization_genes, character()),
  hormone_score = compute_signature_score(expr244_mat, hormone_response_genes, character()),
  lesion_score = compute_signature_score(expr244_mat, lesion_up, lesion_down),
  adeno_consensus_score = compute_signature_score(expr244_mat, adeno_up, adeno_down),
  stringsAsFactors = FALSE
)
proj244 <- project_to_timeline(score244$receptivity_score, calib_model)
all_scores[["adeno_gse244236"]] <- bind_cols(
  data.frame(
    dataset = "GSE244236", disease = "Adenomyosis", cohort_type = "Organoid transcriptomics",
    sample = score244$sample, group = score244$group, tissue_context = score244$tissue_context,
    receptivity_score = score244$receptivity_score,
    timing_state = as.character(proj244$timing_state),
    rif_immune_score = score244$rif_immune_score,
    decidualization_score = score244$decidualization_score,
    hormone_score = score244$hormone_score,
    lesion_score = score244$lesion_score,
    adeno_consensus_score = score244$adeno_consensus_score
  ),
  proj244[, c("pred_day", "woi_distance")]
)

# Adenomyosis tissue GSE190580
expr190 <- as.data.frame(read_excel("/Volumes/Extreme SSD/04_adenomyosis/GSE190580/GSE190580_Raw_count_data.xlsx", sheet = "Gene_COUNTS"))
names(expr190)[1] <- "ensembl"
expr190$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = expr190$ensembl, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
expr190 <- expr190 %>% filter(!is.na(symbol), symbol != "") %>% dplyr::select(-ensembl) %>% group_by(symbol) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
expr190_mat <- as.data.frame(expr190)
rownames(expr190_mat) <- expr190_mat$symbol
expr190_mat$symbol <- NULL
expr190_mat <- log2(as.matrix(expr190_mat) + 1)
colinfo190 <- data.frame(sample = colnames(expr190_mat), stringsAsFactors = FALSE)
colinfo190$group <- ifelse(grepl("^Ade-", colinfo190$sample), "Adenomyosis", "Control")
colinfo190$compartment <- ifelse(grepl("endo", colinfo190$sample, ignore.case = TRUE), "Endometrium", "Myometrium")
score190 <- data.frame(
  sample = colnames(expr190_mat),
  group = colinfo190$group,
  tissue_context = colinfo190$compartment,
  receptivity_score = compute_signature_score(expr190_mat, late_genes, early_genes),
  rif_immune_score = compute_signature_score(expr190_mat, rif_up_genes, rif_down_genes),
  decidualization_score = compute_signature_score(expr190_mat, decidualization_genes, character()),
  hormone_score = compute_signature_score(expr190_mat, hormone_response_genes, character()),
  lesion_score = compute_signature_score(expr190_mat, lesion_up, lesion_down),
  adeno_consensus_score = compute_signature_score(expr190_mat, adeno_up, adeno_down),
  stringsAsFactors = FALSE
)
proj190 <- project_to_timeline(score190$receptivity_score, calib_model)
all_scores[["adeno_gse190580"]] <- bind_cols(
  data.frame(
    dataset = "GSE190580", disease = "Adenomyosis", cohort_type = "Tissue transcriptomics",
    sample = score190$sample, group = score190$group, tissue_context = score190$tissue_context,
    receptivity_score = score190$receptivity_score,
    timing_state = as.character(proj190$timing_state),
    rif_immune_score = score190$rif_immune_score,
    decidualization_score = score190$decidualization_score,
    hormone_score = score190$hormone_score,
    lesion_score = score190$lesion_score,
    adeno_consensus_score = score190$adeno_consensus_score
  ),
  proj190[, c("pred_day", "woi_distance")]
)

# Adenomyosis stromal GSE157718
expr157 <- read.delim(gzfile("/Volumes/Extreme SSD/04_adenomyosis/GSE157718/GSE157718_gene_tpm_matrix.txt.gz"), check.names = FALSE)
names(expr157)[1] <- "ensembl"
expr157$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = expr157$ensembl, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
expr157 <- expr157 %>% filter(!is.na(symbol), symbol != "") %>% dplyr::select(-ensembl) %>% group_by(symbol) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
expr157_mat <- as.data.frame(expr157)
rownames(expr157_mat) <- expr157_mat$symbol
expr157_mat$symbol <- NULL
expr157_mat <- log2(as.matrix(expr157_mat) + 1)
score157 <- data.frame(
  sample = colnames(expr157_mat),
  group = ifelse(grepl("^ES", colnames(expr157_mat)), "Adenomyosis", "Control"),
  tissue_context = "Stroma",
  receptivity_score = compute_signature_score(expr157_mat, late_genes, early_genes),
  rif_immune_score = compute_signature_score(expr157_mat, rif_up_genes, rif_down_genes),
  decidualization_score = compute_signature_score(expr157_mat, decidualization_genes, character()),
  hormone_score = compute_signature_score(expr157_mat, hormone_response_genes, character()),
  lesion_score = compute_signature_score(expr157_mat, lesion_up, lesion_down),
  adeno_consensus_score = compute_signature_score(expr157_mat, adeno_up, adeno_down),
  stringsAsFactors = FALSE
)
proj157 <- project_to_timeline(score157$receptivity_score, calib_model)
all_scores[["adeno_gse157718"]] <- bind_cols(
  data.frame(
    dataset = "GSE157718", disease = "Adenomyosis", cohort_type = "Stromal transcriptomics",
    sample = score157$sample, group = score157$group, tissue_context = score157$tissue_context,
    receptivity_score = score157$receptivity_score,
    timing_state = as.character(proj157$timing_state),
    rif_immune_score = score157$rif_immune_score,
    decidualization_score = score157$decidualization_score,
    hormone_score = score157$hormone_score,
    lesion_score = score157$lesion_score,
    adeno_consensus_score = score157$adeno_consensus_score
  ),
  proj157[, c("pred_day", "woi_distance")]
)

# Adenomyosis whole tissue GSE78851
gse788 <- read_series_matrix("/Volumes/Extreme SSD/04_adenomyosis/GSE78851/GSE78851_series_matrix.txt.gz")
expr788 <- gse788$expr
names(expr788)[1] <- "probe_id"
expr788$symbol <- AnnotationDbi::mapIds(hugene10sttranscriptcluster.db, keys = as.character(expr788$probe_id), keytype = "PROBEID", column = "SYMBOL", multiVals = "first")
expr788 <- expr788 %>% filter(!is.na(symbol), symbol != "") %>% dplyr::select(-probe_id) %>% group_by(symbol) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
expr788_mat <- as.data.frame(expr788)
rownames(expr788_mat) <- expr788_mat$symbol
expr788_mat$symbol <- NULL
score788 <- data.frame(
  sample = colnames(expr788_mat),
  group = ifelse(grepl("Adenomyosis", gse788$meta$title), "Adenomyosis", "Control"),
  tissue_context = "Whole_tissue",
  receptivity_score = compute_signature_score(as.matrix(expr788_mat), late_genes, early_genes),
  rif_immune_score = compute_signature_score(as.matrix(expr788_mat), rif_up_genes, rif_down_genes),
  decidualization_score = compute_signature_score(as.matrix(expr788_mat), decidualization_genes, character()),
  hormone_score = compute_signature_score(as.matrix(expr788_mat), hormone_response_genes, character()),
  lesion_score = compute_signature_score(as.matrix(expr788_mat), lesion_up, lesion_down),
  adeno_consensus_score = compute_signature_score(as.matrix(expr788_mat), adeno_up, adeno_down),
  stringsAsFactors = FALSE
)
proj788 <- project_to_timeline(score788$receptivity_score, calib_model)
all_scores[["adeno_gse78851"]] <- bind_cols(
  data.frame(
    dataset = "GSE78851", disease = "Adenomyosis", cohort_type = "Whole-tissue transcriptomics",
    sample = score788$sample, group = score788$group, tissue_context = score788$tissue_context,
    receptivity_score = score788$receptivity_score,
    timing_state = as.character(proj788$timing_state),
    rif_immune_score = score788$rif_immune_score,
    decidualization_score = score788$decidualization_score,
    hormone_score = score788$hormone_score,
    lesion_score = score788$lesion_score,
    adeno_consensus_score = score788$adeno_consensus_score
  ),
  proj788[, c("pred_day", "woi_distance")]
)

atlas_scores <- bind_rows(all_scores)
write.csv(atlas_scores, file.path(tabdir, "cross_disease_score_table.csv"), row.names = FALSE)

atlas_groups <- atlas_scores %>%
  mutate(group_label = paste(dataset, group, tissue_context, sep = " | ")) %>%
  group_by(disease, dataset, cohort_type, group, tissue_context, group_label) %>%
  summarise(
    n = n(),
    mean_pred_day = mean(pred_day, na.rm = TRUE),
    mean_woi_distance = mean(woi_distance, na.rm = TRUE),
    mean_receptivity = mean(receptivity_score, na.rm = TRUE),
    mean_immune = mean(rif_immune_score, na.rm = TRUE),
    mean_decidualization = mean(decidualization_score, na.rm = TRUE),
    mean_hormone = mean(hormone_score, na.rm = TRUE),
    mean_lesion = mean(lesion_score, na.rm = TRUE),
    mean_adeno_consensus = mean(adeno_consensus_score, na.rm = TRUE),
    prop_delayed = mean(timing_state == "Delayed", na.rm = TRUE),
    prop_in_phase = mean(timing_state == "In-phase", na.rm = TRUE),
    prop_advanced = mean(timing_state == "Advanced", na.rm = TRUE),
    .groups = "drop"
  )
write.csv(atlas_groups, file.path(tabdir, "cross_disease_group_summary.csv"), row.names = FALSE)

# Figure 1: projected day by disease
plot_groups <- atlas_scores %>%
  filter(
    (disease == "RIF" & group %in% c("Fertile_LH7", "RIF", "Control")) |
      (disease == "Endometriosis" & group %in% c("Control_Eutopic", "Endo_Eutopic", "Ectopic", "endometriosis", "Control", "patient with endometriosis", "eutopic endometrium", "endometriotic lesion")) |
      (disease == "Adenomyosis" & group %in% c("Control", "Adenomyosis"))
  ) %>%
  mutate(plot_group = case_when(
    group %in% c("Fertile_LH7", "Control", "Control_Eutopic", "eutopic endometrium") ~ "Control-like",
    group %in% c("Endo_Eutopic", "endometriosis", "patient with endometriosis") ~ "Disease eutopic",
    group %in% c("Ectopic", "endometriotic lesion") ~ "Lesion",
    TRUE ~ group
  ))

p1 <- ggplot(plot_groups, aes(plot_group, pred_day, color = plot_group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1.1) +
  facet_wrap(~ disease, scales = "free_x", nrow = 1) +
  geom_hline(yintercept = 7, linetype = 2, color = "grey45") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = "Cross-disease projection onto the fertile receptivity timeline", x = "", y = "Predicted receptive day")
ggsave(file.path(figdir, "Figure_ATLAS_1_cross_disease_projection.png"), p1, width = 12.8, height = 5.8, dpi = 260)

# Figure 2: timing state fractions
state_df <- atlas_scores %>%
  filter(!is.na(timing_state)) %>%
  mutate(group_simple = case_when(
    disease == "RIF" & group %in% c("RIF", "Control", "Fertile_LH7") ~ group,
    disease == "Endometriosis" & group %in% c("Control_Eutopic", "Endo_Eutopic", "Ectopic", "endometriosis", "Control", "patient with endometriosis", "eutopic endometrium", "endometriotic lesion") ~ group,
    disease == "Adenomyosis" ~ paste(group, tissue_context, sep = "_"),
    TRUE ~ paste(group, tissue_context, sep = "_")
  )) %>%
  mutate(group_label = case_when(
    group_simple == "Adenomyosis_Endometrium" ~ "Adenomyosis endometrium",
    group_simple == "Adenomyosis_GEST" ~ "Adenomyosis GEST",
    group_simple == "Adenomyosis_Myometrium" ~ "Adenomyosis myometrium",
    group_simple == "Adenomyosis_SEC" ~ "Adenomyosis SEC",
    group_simple == "Adenomyosis_Stroma" ~ "Adenomyosis stroma",
    group_simple == "Adenomyosis_Whole_tissue" ~ "Adenomyosis whole tissue",
    group_simple == "Control_Endometrium" ~ "Control endometrium",
    group_simple == "Control_GEST" ~ "Control GEST",
    group_simple == "Control_Myometrium" ~ "Control myometrium",
    group_simple == "Control_SEC" ~ "Control SEC",
    group_simple == "Control_Stroma" ~ "Control stroma",
    group_simple == "Control_Whole_tissue" ~ "Control whole tissue",
    group_simple == "Control_Eutopic" ~ "Control eutopic",
    group_simple == "Endo_Eutopic" ~ "Endo eutopic",
    group_simple == "Ectopic" ~ "Ectopic lesion",
    group_simple == "Ectopic_Adjacent_Ectopic" ~ "Adjacent ectopic",
    group_simple == "Ectopic_Ovary_Ectopic_Ovary" ~ "Ovarian ectopic",
    group_simple == "endometriosis" ~ "Bulk endometriosis",
    group_simple == "endometriotic lesion" ~ "Lesion",
    group_simple == "eutopic endometrium" ~ "Eutopic endometrium",
    group_simple == "healthy control_linearly control" ~ "Healthy control",
    group_simple == "patient with endometriosis" ~ "Patient endometriosis",
    group_simple == "Control" ~ "Control",
    group_simple == "Fertile_LH7" ~ "Fertile LH7",
    group_simple == "LH11_Fertile_reference" ~ "LH11 fertile ref",
    group_simple == "LH3_Fertile_reference" ~ "LH3 fertile ref",
    group_simple == "LH5_Fertile_reference" ~ "LH5 fertile ref",
    group_simple == "LH9_Fertile_reference" ~ "LH9 fertile ref",
    group_simple == "RIF" ~ "RIF",
    TRUE ~ gsub("_", " ", group_simple)
  )) %>%
  group_by(disease, group_simple, group_label, timing_state) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  group_by(disease, group_simple, group_label) %>%
  mutate(group_label = factor(group_label, levels = unique(group_label))) %>%
  ungroup()
p2 <- ggplot(state_df, aes(group_label, prop, fill = timing_state)) +
  geom_col() +
  facet_wrap(~ disease, scales = "free_x", nrow = 1) +
  theme_bw(base_size = 12) +
  coord_cartesian(clip = "off") +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8.2, lineheight = 0.95),
    plot.margin = margin(10, 28, 96, 30),
    panel.spacing.x = grid::unit(10, "pt")
  ) +
  labs(title = "Timing-state decomposition across implantation-failure-associated disorders", x = "", y = "Proportion")
ggsave(file.path(figdir, "Figure_ATLAS_2_timing_state_fraction.png"), p2, width = 17.8, height = 7.8, dpi = 260)

# Figure 3: module heatmap
heat_groups <- atlas_groups %>%
  filter(
    (disease == "RIF" & group %in% c("Fertile_LH7", "RIF")) |
      (disease == "Endometriosis" & group %in% c("Control_Eutopic", "Endo_Eutopic", "Ectopic", "Control", "endometriosis", "patient with endometriosis", "eutopic endometrium", "endometriotic lesion")) |
      (disease == "Adenomyosis" & group %in% c("Control", "Adenomyosis"))
  ) %>%
  mutate(label = paste(disease, dataset, group, tissue_context, sep = " | ")) %>%
  dplyr::select(
    label,
    mean_woi_distance,
    mean_immune,
    mean_decidualization,
    mean_hormone,
    mean_lesion,
    mean_adeno_consensus
  )
heat_mat <- as.matrix(heat_groups[, -1])
rownames(heat_mat) <- heat_groups$label
heat_mat <- scale(heat_mat)
png(file.path(figdir, "Figure_ATLAS_3_cross_disease_module_heatmap.png"), width = 2400, height = 1800, res = 240)
pheatmap(t(heat_mat), cluster_cols = TRUE, cluster_rows = FALSE, main = "Cross-disease molecular module atlas")
dev.off()

# Figure 4: spatial RIF niche
rif_spatial <- atlas_scores %>% filter(dataset == "GSE287278")
cor_rif <- suppressWarnings(cor(rif_spatial$receptivity_score, rif_spatial$rif_immune_score, method = "spearman", use = "pairwise.complete.obs"))
p4 <- ggplot(rif_spatial, aes(receptivity_score, rif_immune_score, color = group)) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  theme_bw(base_size = 12) +
  labs(
    title = paste0("RIF spatial niche organization (Spearman rho = ", round(cor_rif, 3), ")"),
    x = "Receptivity score",
    y = "Immune-activation score"
  )
ggsave(file.path(figdir, "Figure_ATLAS_4_rif_spatial_niche.png"), p4, width = 7.8, height = 5.8, dpi = 260)

# Figure 5: endometriosis spatial niche
endo_spatial <- atlas_scores %>% filter(dataset == "GSE263897")
endo_spatial$niche <- ifelse(
  endo_spatial$receptivity_score < median(endo_spatial$receptivity_score, na.rm = TRUE) &
    endo_spatial$rif_immune_score > median(endo_spatial$rif_immune_score, na.rm = TRUE),
  "Immune-high / receptivity-low", "Other"
)
p5a <- ggplot(endo_spatial, aes(tissue_context, receptivity_score, color = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75), size = 1.2, alpha = 0.7) +
  theme_bw(base_size = 12) +
  labs(title = "Endometriosis spatial receptivity by ROI cell type", x = "", y = "Receptivity score")
ggsave(file.path(figdir, "Figure_ATLAS_5A_endo_spatial_receptivity.png"), p5a, width = 9.2, height = 5.6, dpi = 260)
p5b <- ggplot(endo_spatial, aes(tissue_context, rif_immune_score, color = group)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75), size = 1.2, alpha = 0.7) +
  theme_bw(base_size = 12) +
  labs(title = "Endometriosis spatial immune-activation by ROI cell type", x = "", y = "Immune-activation score")
ggsave(file.path(figdir, "Figure_ATLAS_5B_endo_spatial_immune.png"), p5b, width = 9.2, height = 5.6, dpi = 260)

# Figure 6: cross-dataset effect sizes for WOI distance
effect_specs <- list(
  list(dataset = "GSE111974", disease = "RIF", case = c("RIF"), ctrl = c("Control")),
  list(dataset = "GSE58144", disease = "RIF", case = c("RIF"), ctrl = c("Control")),
  list(dataset = "GSE179640", disease = "Endometriosis", case = c("Endo_Eutopic"), ctrl = c("Control_Eutopic")),
  list(dataset = "GSE214411", disease = "Endometriosis", case = c("endometriosis"), ctrl = c("Control")),
  list(dataset = "GSE135485", disease = "Endometriosis", case = c("patient with endometriosis"), ctrl = c("healthy control")),
  list(dataset = "GSE244236", disease = "Adenomyosis", case = c("Adenomyosis"), ctrl = c("Control")),
  list(dataset = "GSE190580", disease = "Adenomyosis", case = c("Adenomyosis"), ctrl = c("Control")),
  list(dataset = "GSE157718", disease = "Adenomyosis", case = c("Adenomyosis"), ctrl = c("Control")),
  list(dataset = "GSE78851", disease = "Adenomyosis", case = c("Adenomyosis"), ctrl = c("Control"))
)

effect_df <- bind_rows(lapply(effect_specs, function(spec) {
  dat <- atlas_scores %>% filter(dataset == spec$dataset, group %in% c(spec$case, spec$ctrl))
  if (spec$dataset == "GSE190580") dat <- dat %>% filter(tissue_context == "Endometrium")
  grp <- ifelse(dat$group %in% spec$case, 1, 0)
  eff <- compute_effect_size(dat$woi_distance, grp)
  cbind(data.frame(dataset = spec$dataset, disease = spec$disease), eff)
}))
write.csv(effect_df, file.path(tabdir, "cross_disease_woi_effect_sizes.csv"), row.names = FALSE)
meta_df <- effect_df %>% group_by(disease) %>% group_modify(~ random_effects_meta(.x)) %>% ungroup()
write.csv(meta_df, file.path(tabdir, "cross_disease_random_effects_summary.csv"), row.names = FALSE)
leave_one_out <- bind_rows(lapply(seq_len(nrow(effect_df)), function(i) {
  tmp <- effect_df[-i, , drop = FALSE]
  out <- tmp %>% group_by(disease) %>% group_modify(~ random_effects_meta(.x)) %>% ungroup()
  out$left_out <- effect_df$dataset[i]
  out
}))
write.csv(leave_one_out, file.path(tabdir, "cross_disease_leave_one_out_meta.csv"), row.names = FALSE)

p6 <- ggplot(effect_df, aes(effect, reorder(paste(disease, dataset, sep = " | "), effect), color = disease)) +
  geom_point(size = 2.4) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.18, orientation = "y") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
  theme_bw(base_size = 12) +
  labs(title = "Cross-dataset effect sizes for window-of-implantation displacement", x = "Standardized mean difference of WOI distance", y = "")
ggsave(file.path(figdir, "Figure_ATLAS_6_woi_effect_forest.png"), p6, width = 10.2, height = 6.8, dpi = 260)

summary_df <- data.frame(
  n_total_rows = nrow(atlas_scores),
  n_datasets = length(unique(atlas_scores$dataset)),
  n_diseases = length(unique(atlas_scores$disease)),
  rif_mean_pred_day = mean(atlas_scores$pred_day[atlas_scores$disease == "RIF" & atlas_scores$group == "RIF"], na.rm = TRUE),
  endo_eutopic_mean_pred_day = mean(atlas_scores$pred_day[atlas_scores$dataset == "GSE179640" & atlas_scores$group == "Endo_Eutopic"], na.rm = TRUE),
  adeno_organoid_mean_pred_day = mean(atlas_scores$pred_day[atlas_scores$dataset == "GSE244236" & atlas_scores$group == "Adenomyosis"], na.rm = TRUE),
  rif_spatial_correlation = cor_rif
)
write.csv(summary_df, file.path(tabdir, "implantation_failure_atlas_summary.csv"), row.names = FALSE)

print(summary_df)