#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(readxl)
  library(Matrix)
  library(limma)
  library(fgsea)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

set.seed(123)
Sys.setenv(HOME = "/private/tmp", XDG_CACHE_HOME = "/private/tmp")

root <- "[local path omitted]"
outdir <- file.path(root, "analysis/05_implantation_failure_atlas/reviewer_stats")
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")
msigdb_zip <- "[local path omitted]"
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

project_day_states <- function(pred_day) {
  cut(pred_day, breaks = c(-Inf, 6.5, 8.5, Inf), labels = c("Delayed", "In-phase", "Advanced"))
}

scale_in_dataset <- function(x) {
  z <- as.numeric(scale(x))
  z[!is.finite(z)] <- NA_real_
  z
}

bootstrap_mean_ci <- function(x, nboot = 3000) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(mean = mean(x), lower = NA_real_, upper = NA_real_))
  boots <- replicate(nboot, mean(sample(x, length(x), replace = TRUE)))
  c(mean = mean(x), lower = unname(quantile(boots, 0.025)), upper = unname(quantile(boots, 0.975)))
}

wrap_text <- function(x, width = 28) {
  vapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n", character(1))
}

read_gmt_simple <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- vector("list", length(lines))
  nms <- character(length(lines))
  for (i in seq_along(lines)) {
    parts <- strsplit(lines[[i]], "\t", fixed = TRUE)[[1]]
    nms[[i]] <- parts[[1]]
    out[[i]] <- unique(parts[-c(1, 2)])
  }
  names(out) <- nms
  out
}

compute_effect_size <- function(score, group01) {
  idx1 <- which(group01 == 1)
  idx0 <- which(group01 == 0)
  s1 <- score[idx1]
  s0 <- score[idx0]
  m1 <- mean(s1, na.rm = TRUE)
  m0 <- mean(s0, na.rm = TRUE)
  sd1 <- sd(s1, na.rm = TRUE)
  sd0 <- sd(s0, na.rm = TRUE)
  n1 <- length(idx1)
  n0 <- length(idx0)
  sp <- sqrt(((n1 - 1) * sd1^2 + (n0 - 1) * sd0^2) / pmax(n1 + n0 - 2, 1))
  d <- (m1 - m0) / sp
  se <- sqrt((n1 + n0) / (n1 * n0) + (d^2 / (2 * pmax(n1 + n0 - 2, 1))))
  data.frame(
    effect = d,
    se = se,
    lower = d - 1.96 * se,
    upper = d + 1.96 * se,
    n_case = n1,
    n_ctrl = n0
  )
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
  data.frame(
    effect = pooled,
    se = se_pooled,
    lower = pooled - 1.96 * se_pooled,
    upper = pooled + 1.96 * se_pooled,
    tau2 = tau2
  )
}

exact_signflip_p <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n == 0) return(NA_real_)
  if (n > 15) return(NA_real_)
  obs <- abs(mean(x))
  signs <- expand.grid(rep(list(c(-1, 1)), n))
  vals <- apply(signs, 1, function(s) abs(mean(as.numeric(s) * abs(x))))
  mean(vals >= obs)
}

run_gsego_table <- function(stats, label) {
  stats <- stats[is.finite(stats)]
  eg <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = names(stats), keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
  keep <- !is.na(eg)
  stats2 <- stats[keep]
  names(stats2) <- unname(eg[keep])
  stats2 <- sort(stats2, decreasing = TRUE)
  stats2 <- stats2[!duplicated(names(stats2))]
  gsea <- suppressWarnings(
    clusterProfiler::gseGO(
      geneList = stats2,
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      keyType = "ENTREZID",
      minGSSize = 15,
      maxGSSize = 500,
      pvalueCutoff = 1,
      verbose = FALSE
    )
  )
  fg <- as.data.frame(gsea@result)
  if (nrow(fg) == 0) {
    return(data.frame())
  }
  fg %>%
    arrange(p.adjust, desc(abs(NES))) %>%
    mutate(contrast = label, collection = "GO_BP")
}

run_fgsea_gmt <- function(stats, pathways, label, collection_label = "Hallmark") {
  stats <- stats[is.finite(stats)]
  stats <- sort(stats, decreasing = TRUE)
  fg <- suppressWarnings(
    fgsea::fgsea(
      pathways = pathways,
      stats = stats,
      minSize = 10,
      maxSize = 500
    )
  )
  fg <- as.data.frame(fg)
  if (nrow(fg) == 0) {
    return(data.frame())
  }
  fg %>%
    arrange(padj, desc(abs(NES))) %>%
    mutate(
      contrast = label,
      collection = collection_label,
      p.adjust = padj,
      Description = pathway
    )
}

read_series_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  table_start <- which(lines == "!series_matrix_table_begin")
  table_end <- which(lines == "!series_matrix_table_end")
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

atlas <- read.csv(file.path(root, "analysis/05_implantation_failure_atlas/final/tables/cross_disease_score_table.csv"), check.names = FALSE)
atlas_ext <- atlas %>%
  group_by(dataset) %>%
  mutate(
    receptivity_low_z = scale_in_dataset(-receptivity_score),
    immune_high_z = scale_in_dataset(rif_immune_score),
    decidual_low_z = scale_in_dataset(-decidualization_score),
    implantation_module_score = rowMeans(cbind(receptivity_low_z, immune_high_z, decidual_low_z), na.rm = TRUE),
    common_component_n = rowSums(!is.na(cbind(receptivity_low_z, immune_high_z, decidual_low_z)))
  ) %>%
  ungroup() %>%
  mutate(
    implantation_module_score = ifelse(common_component_n >= 2, implantation_module_score, NA_real_),
    ablation_drop_receptivity = ifelse(rowSums(!is.na(cbind(immune_high_z, decidual_low_z))) == 2, rowMeans(cbind(immune_high_z, decidual_low_z), na.rm = TRUE), NA_real_),
    ablation_drop_immune = ifelse(rowSums(!is.na(cbind(receptivity_low_z, decidual_low_z))) == 2, rowMeans(cbind(receptivity_low_z, decidual_low_z), na.rm = TRUE), NA_real_),
    ablation_drop_decidual = ifelse(rowSums(!is.na(cbind(receptivity_low_z, immune_high_z))) == 2, rowMeans(cbind(receptivity_low_z, immune_high_z), na.rm = TRUE), NA_real_)
  )
write.csv(atlas_ext, file.path(tabdir, "cross_disease_score_table_with_common_axis.csv"), row.names = FALSE)

plot_unit_map <- c(
  "Single-cell pseudobulk" = "scRNA pseudobulk sample",
  "Bulk endometrium" = "bulk sample",
  "Tissue transcriptomics" = "bulk sample",
  "Whole-tissue transcriptomics" = "bulk sample",
  "Stromal transcriptomics" = "bulk sample",
  "Compartment-aware tissue transcriptomics" = "bulk sample",
  "Organoid transcriptomics" = "organoid sample",
  "Spatial transcriptomics" = "Visium section",
  "GeoMx spatial transcriptomics" = "GeoMx ROI"
)

figure2_tbl <- atlas_ext %>%
  mutate(
    plot_unit = recode(cohort_type, !!!plot_unit_map, .default = cohort_type),
    plot_unit = ifelse(is.na(plot_unit) | plot_unit == "NA", "bulk sample", plot_unit),
    section_id = ifelse(dataset == "GSE287278", sub("::.*$", "", sample), sample),
    plot_id = ifelse(dataset == "GSE287278", section_id, sample)
  ) %>%
  filter(!is.na(pred_day)) %>%
  group_by(dataset, disease, cohort_type, plot_unit, group, tissue_context, plot_id) %>%
  summarise(pred_day = mean(pred_day, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    plot_unit = factor(plot_unit, levels = c("bulk sample", "scRNA pseudobulk sample", "Visium section", "GeoMx ROI", "organoid sample")),
    disease = factor(disease, levels = c("RIF", "Endometriosis", "Adenomyosis"))
  )

p_units <- ggplot(figure2_tbl, aes(disease, pred_day, color = disease)) +
  geom_boxplot(outlier.shape = NA, width = 0.52, alpha = 0.18) +
  geom_jitter(width = 0.14, alpha = 0.55, size = 1.7) +
  facet_wrap(~ plot_unit, scales = "free_y", ncol = 2) +
  theme_bw(base_size = 12) +
  guides(color = "none") +
  labs(
    title = "Cross-disease projection across observation units",
    subtitle = "Points represent sample-, section-, ROI-, or organoid-level analytical observations rather than all independent patients",
    x = "",
    y = "Predicted receptive day"
  )
ggsave(file.path(figdir, "Figure_ATLAS_RS_14_cross_disease_projection_units.png"), p_units, width = 9.8, height = 7.4, dpi = 260)

module_specs <- c(
  "woi_distance",
  "rif_immune_score",
  "lesion_score",
  "adeno_consensus_score",
  "hormone_score"
)

case_only <- atlas_ext %>%
  filter(
    (disease == "RIF" & group == "RIF") |
      (disease == "Endometriosis" & !(group %in% c("Control_Eutopic", "Control", "healthy control", "No_endo_detected"))) |
      (disease == "Adenomyosis" & group == "Adenomyosis")
  )

atlas_ext <- atlas_ext %>%
  mutate(case_flag = case_when(
    dataset == "GSE250130" ~ ifelse(group == "RIF", 1, ifelse(group == "Fertile_LH7", 0, NA_real_)),
    dataset %in% c("GSE111974", "GSE58144", "GSE287278") ~ ifelse(group == "RIF", 1, ifelse(group == "Control", 0, NA_real_)),
    dataset == "GSE179640" ~ ifelse(group %in% c("Endo_Eutopic", "Ectopic", "Ectopic_Adjacent", "Ectopic_Ovary", "Mix"), 1, ifelse(group == "Control_Eutopic", 0, NA_real_)),
    dataset == "GSE213216" ~ ifelse(group %in% c("Lesion", "Eutopic", "Ovary"), 1, ifelse(group == "No_endo_detected", 0, NA_real_)),
    dataset == "GSE214411" ~ ifelse(group == "endometriosis", 1, ifelse(group == "Control", 0, NA_real_)),
    dataset == "GSE135485" ~ ifelse(group == "patient with endometriosis", 1, ifelse(group == "healthy control", 0, NA_real_)),
    dataset == "GSE263897" ~ ifelse(group == "endometriotic lesion", 1, ifelse(group == "eutopic endometrium", 0, NA_real_)),
    dataset %in% c("GSE244236", "GSE190580", "GSE157718", "GSE78851") ~ ifelse(group == "Adenomyosis", 1, ifelse(group == "Control", 0, NA_real_)),
    TRUE ~ NA_real_
  ))

derive_dataset_effects <- function(input_df, score_col = "implantation_module_score") {
  bind_rows(lapply(module_specs, function(mod) {
  bind_rows(lapply(sort(unique(na.omit(input_df$dataset))), function(ds) {
    dat <- input_df %>%
      filter(dataset == ds, !is.na(case_flag), !is.na(.data[[score_col]]), !is.na(.data[[mod]]))
    if (nrow(dat) < 6 || length(unique(dat$case_flag)) < 2) return(NULL)
    fit <- lm(as.formula(paste0(mod, " ~ ", score_col)), data = dat)
    dat$residual_value <- resid(fit)
    eff <- compute_effect_size(dat$residual_value, dat$case_flag)
    data.frame(
      disease = unique(dat$disease),
      dataset = ds,
      module = mod,
      effect = eff$effect,
      se = eff$se,
      lower = eff$lower,
      upper = eff$upper,
      n_case = eff$n_case,
      n_ctrl = eff$n_ctrl
    )
  }))
  }))
}

dataset_effects <- derive_dataset_effects(atlas_ext)
write.csv(dataset_effects, file.path(tabdir, "disease_specific_residual_dataset_effects.csv"), row.names = FALSE)

residual_summary <- dataset_effects %>%
  group_by(disease, module) %>%
  group_modify(~ {
    pooled <- random_effects_meta(.x)
    data.frame(
      effect = pooled$effect,
      se = pooled$se,
      lower = pooled$lower,
      upper = pooled$upper,
      tau2 = pooled$tau2,
      n_datasets = nrow(.x)
    )
  }) %>%
  ungroup() %>%
  mutate(
    module_label = recode(
      module,
      woi_distance = "Timing instability",
      rif_immune_score = "Immune residual",
      lesion_score = "Lesion program",
      adeno_consensus_score = "Adenomyosis program",
      hormone_score = "Hormone-response disruption"
    )
  )
write.csv(residual_summary, file.path(tabdir, "disease_specific_residual_module_summary.csv"), row.names = FALSE)

atlas_complete <- atlas_ext %>% filter(common_component_n == 3)
dataset_effects_complete <- derive_dataset_effects(atlas_complete)
write.csv(dataset_effects_complete, file.path(tabdir, "common_axis_complete_component_dataset_effects.csv"), row.names = FALSE)

residual_summary_complete <- dataset_effects_complete %>%
  group_by(disease, module) %>%
  group_modify(~ {
    pooled <- random_effects_meta(.x)
    data.frame(
      effect_complete = pooled$effect,
      se_complete = pooled$se,
      lower_complete = pooled$lower,
      upper_complete = pooled$upper,
      tau2_complete = pooled$tau2,
      n_datasets_complete = nrow(.x)
    )
  }) %>%
  ungroup()

complete_compare <- residual_summary %>%
  dplyr::select(disease, module, module_label, effect_primary = effect, lower_primary = lower, upper_primary = upper, n_datasets_primary = n_datasets) %>%
  left_join(residual_summary_complete, by = c("disease", "module")) %>%
  mutate(direction_consistent = sign(effect_primary) == sign(effect_complete))
write.csv(complete_compare, file.path(tabdir, "common_axis_complete_component_sensitivity.csv"), row.names = FALSE)

complete_plot <- complete_compare %>%
  mutate(
    module_label = factor(module_label, levels = c("Timing instability", "Immune residual", "Lesion program", "Adenomyosis program", "Hormone-response disruption")),
    disease = factor(disease, levels = c("RIF", "Endometriosis", "Adenomyosis"))
  ) %>%
  dplyr::select(disease, module_label, effect_primary, effect_complete) %>%
  pivot_longer(cols = c(effect_primary, effect_complete), names_to = "analysis", values_to = "effect") %>%
  mutate(analysis = recode(analysis, effect_primary = "Primary (>=2 components)", effect_complete = "Sensitivity (3 components only)"))

p_complete <- ggplot(complete_plot, aes(module_label, disease, fill = effect)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", effect)), size = 3.1) +
  facet_wrap(~ analysis) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1)) +
  labs(
    title = "Common-axis component sensitivity of disease-specific residual programs",
    x = "",
    y = "",
    fill = "Pooled\nresidual effect"
  )
ggsave(file.path(figdir, "Figure_ATLAS_RS_15_complete_component_sensitivity.png"), p_complete, width = 10.2, height = 4.8, dpi = 260)

p_common <- case_only %>%
  filter(!is.na(implantation_module_score)) %>%
  ggplot(aes(disease, implantation_module_score, color = disease)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.14, alpha = 0.45, size = 1.5) +
  theme_bw(base_size = 12) +
  labs(
    title = "Hypothesis-driven implantation-failure module score across disorders",
    x = "",
    y = "Implantation-failure module score"
  ) +
  guides(color = "none")
ggsave(file.path(figdir, "Figure_ATLAS_RS_8_common_failure_axis.png"), p_common, width = 7.2, height = 5.2, dpi = 260)

p_resid <- residual_summary %>%
  mutate(
    disease = factor(disease, levels = c("RIF", "Endometriosis", "Adenomyosis")),
    module_label = factor(module_label, levels = c("Timing instability", "Immune residual", "Lesion program", "Adenomyosis program", "Hormone-response disruption"))
  ) %>%
  ggplot(aes(module_label, disease, fill = effect)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", effect)), size = 3.2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 22, hjust = 1)) +
  labs(
    title = "Disease-specific residual programs after removing the shared failure axis",
    x = "",
    y = "",
    fill = "Residual\npooled effect"
  )
ggsave(file.path(figdir, "Figure_ATLAS_RS_9_disease_specific_residual_heatmap.png"), p_resid, width = 8.6, height = 4.6, dpi = 260)

ablation_specs <- c(
  "Primary (3 modules)" = "implantation_module_score",
  "Drop receptivity-low" = "ablation_drop_receptivity",
  "Drop immune-high" = "ablation_drop_immune",
  "Drop decidualization-low" = "ablation_drop_decidual"
)

ablation_summary <- bind_rows(lapply(names(ablation_specs), function(lbl) {
  score_col <- ablation_specs[[lbl]]
  dat_eff <- derive_dataset_effects(atlas_ext, score_col = score_col)
  if (is.null(dat_eff) || nrow(dat_eff) == 0) return(NULL)
  dat_eff %>%
    group_by(disease, module) %>%
    group_modify(~ {
      pooled <- random_effects_meta(.x)
      data.frame(
        effect = pooled$effect,
        lower = pooled$lower,
        upper = pooled$upper,
        n_datasets = nrow(.x)
      )
    }) %>%
    ungroup() %>%
    mutate(
      ablation = lbl,
      module_label = recode(
        module,
        woi_distance = "Timing instability",
        rif_immune_score = "Immune residual",
        lesion_score = "Lesion program",
        adeno_consensus_score = "Adenomyosis program",
        hormone_score = "Hormone-response disruption"
      )
    )
}))

ablation_summary <- ablation_summary %>%
  group_by(disease, module) %>%
  mutate(
    primary_effect = effect[ablation == "Primary (3 modules)"][1],
    direction_consistent = sign(effect) == sign(primary_effect)
  ) %>%
  ungroup()
write.csv(ablation_summary, file.path(tabdir, "implantation_failure_module_ablation_summary.csv"), row.names = FALSE)

p_ablation <- ablation_summary %>%
  mutate(
    disease = factor(disease, levels = c("RIF", "Endometriosis", "Adenomyosis")),
    module_label = factor(module_label, levels = c("Timing instability", "Immune residual", "Lesion program", "Adenomyosis program", "Hormone-response disruption")),
    ablation = factor(ablation, levels = names(ablation_specs))
  ) %>%
  ggplot(aes(module_label, disease, fill = effect)) +
  geom_tile(color = "white") +
  geom_point(aes(shape = direction_consistent), size = 2.5, color = "black", stroke = 0.7) +
  geom_text(aes(label = sprintf("%.2f", effect)), size = 2.9, vjust = -0.95) +
  facet_wrap(~ ablation, ncol = 2) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), labels = c("FALSE" = "Direction changed", "TRUE" = "Direction preserved")) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 24, hjust = 1), legend.position = "bottom") +
  labs(
    title = "Leave-one-module-out sensitivity of disease-specific residual programs",
    subtitle = "Point fill shows pooled residual effect; point shape shows whether direction is preserved relative to the primary three-module score",
    x = "",
    y = "",
    fill = "Residual\npooled effect",
    shape = ""
  )
ggsave(file.path(figdir, "Figure_ATLAS_RS_22_module_ablation.png"), p_ablation, width = 11.2, height = 7.2, dpi = 260)

receptivity_sig <- read.csv(file.path(root, "analysis/03_rif/round5_maximal/tables/receptivity_signature_genes_round5.csv"), check.names = FALSE)
late_genes <- na.omit(receptivity_sig$late_genes)
early_genes <- na.omit(receptivity_sig$early_genes)
rif_sig <- read.csv(file.path(root, "analysis/03_rif/round5_maximal/tables/rif_directional_signature_genes_round5.csv"), check.names = FALSE)
rif_up <- na.omit(rif_sig$rif_up)
rif_down <- na.omit(rif_sig$rif_down)
adeno_sig <- read.csv(file.path(root, "analysis/04_adenomyosis/round4_deep/tables/consensus_adenomyosis_signature.csv"), check.names = FALSE)
adeno_up <- na.omit(adeno_sig$consensus_up)
adeno_down <- na.omit(adeno_sig$consensus_down)
decidual_genes <- c("IGFBP1", "PRL", "LEFTY2", "FOXO1", "WNT4", "HAND2", "IL15", "SPP1", "PAEP", "GPX3")
lesion_deg <- read.csv(file.path(root, "analysis/02_endometriosis/round4_deep/tables/Figure_ENDO_2A_Ectopic_vs_Eutopic_deg.csv"), check.names = FALSE)
lesion_up <- lesion_deg %>% filter(adj.P.Val < 0.05, logFC > 0) %>% arrange(adj.P.Val) %>% slice_head(n = 30) %>% pull(gene)
lesion_down <- lesion_deg %>% filter(adj.P.Val < 0.05, logFC < 0) %>% arrange(adj.P.Val) %>% slice_head(n = 30) %>% pull(gene)

gene_level_residual_work <- list()
hallmark_work <- list()

hallmark_gmt <- NULL
if (file.exists(msigdb_zip)) {
  zip_members <- utils::unzip(msigdb_zip, list = TRUE)
  hallmark_member <- zip_members$Name[grepl("h\\.all\\.v2026\\.1\\.Hs\\.symbols\\.gmt$", zip_members$Name)]
  if (length(hallmark_member) > 0) {
    hallmark_tmp_dir <- tempfile("hallmark_gmt_")
    dir.create(hallmark_tmp_dir)
    utils::unzip(msigdb_zip, files = hallmark_member[1], exdir = hallmark_tmp_dir)
    hallmark_path <- file.path(hallmark_tmp_dir, hallmark_member[1])
    hallmark_gmt <- read_gmt_simple(hallmark_path)
  }
}

# RIF representative dataset: GSE250130
expr_rif <- read.csv(file.path(root, "analysis/03_rif/tables/GSE250130_firstpass/pseudobulk_counts.csv"), check.names = FALSE)
rownames(expr_rif) <- expr_rif$gene
expr_rif$gene <- NULL
expr_rif <- log2(as.matrix(expr_rif) + 1)
meta_rif <- atlas %>% filter(dataset == "GSE250130") %>% dplyr::select(sample, group)
meta_rif <- meta_rif[match(colnames(expr_rif), meta_rif$sample), ]
common_rif <- rowMeans(cbind(
  scale_in_dataset(-compute_signature_score(expr_rif, late_genes, early_genes)),
  scale_in_dataset(compute_signature_score(expr_rif, rif_up, rif_down)),
  scale_in_dataset(-compute_signature_score(expr_rif, decidual_genes, character()))
), na.rm = TRUE)
design_rif <- model.matrix(~ common_rif + I(meta_rif$group == "RIF"))
colnames(design_rif)[3] <- "case"
fit_rif <- eBayes(lmFit(expr_rif, design_rif))
tt_rif <- topTable(fit_rif, coef = "case", number = Inf, sort.by = "none")
stats_rif <- tt_rif$t
names(stats_rif) <- rownames(tt_rif)
gene_level_residual_work[["RIF_GO"]] <- run_gsego_table(stats_rif, "RIF residual deviation")
if (!is.null(hallmark_gmt)) hallmark_work[["RIF_H"]] <- run_fgsea_gmt(stats_rif, hallmark_gmt, "RIF residual deviation")

# Endometriosis representative dataset: GSE179640
expr_endo <- read.csv(file.path(root, "analysis/02_endometriosis/tables/GSE179640_firstpass/pseudobulk_counts.csv"), check.names = FALSE)
rownames(expr_endo) <- expr_endo$gene
expr_endo$gene <- NULL
expr_endo <- log2(as.matrix(expr_endo) + 1)
meta_endo <- read.csv(file.path(root, "analysis/02_endometriosis/metadata/GSE179640_metadata_round2.csv"), check.names = FALSE)
meta_endo <- data.frame(sample = colnames(expr_endo), group = meta_endo$analysis_group[match(colnames(expr_endo), meta_endo$sample)], stringsAsFactors = FALSE)
common_endo <- rowMeans(cbind(
  scale_in_dataset(-compute_signature_score(expr_endo, late_genes, early_genes)),
  scale_in_dataset(compute_signature_score(expr_endo, rif_up, rif_down)),
  scale_in_dataset(-compute_signature_score(expr_endo, decidual_genes, character()))
), na.rm = TRUE)
case_endo <- !(meta_endo$group %in% c("Control_Eutopic"))
design_endo <- model.matrix(~ common_endo + case_endo)
colnames(design_endo)[3] <- "case"
fit_endo <- eBayes(lmFit(expr_endo, design_endo))
tt_endo <- topTable(fit_endo, coef = "case", number = Inf, sort.by = "none")
stats_endo <- tt_endo$t
names(stats_endo) <- rownames(tt_endo)
gene_level_residual_work[["Endometriosis_GO"]] <- run_gsego_table(stats_endo, "Endometriosis residual deviation")
if (!is.null(hallmark_gmt)) hallmark_work[["Endometriosis_H"]] <- run_fgsea_gmt(stats_endo, hallmark_gmt, "Endometriosis residual deviation")

# Adenomyosis representative dataset: GSE244236
expr244_raw <- as.data.frame(read_excel("/Volumes/Extreme SSD/04_adenomyosis/GSE244236/GSE244236_Normalized_counts.xlsx"))
gene_ids244 <- as.character(expr244_raw[[1]])
expr244 <- expr244_raw[, -1, drop = FALSE]
expr244[] <- lapply(expr244, as.numeric)
rownames(expr244) <- gene_ids244
id_map244 <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = rownames(expr244), keytype = "ENTREZID", column = "SYMBOL", multiVals = "first")
expr244$symbol <- unname(id_map244)
expr244 <- expr244 %>% filter(!is.na(symbol), symbol != "")
expr244 <- expr244 %>% dplyr::select(symbol, everything()) %>% group_by(symbol) %>% summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
expr244_mat <- as.data.frame(expr244)
rownames(expr244_mat) <- expr244_mat$symbol
expr244_mat$symbol <- NULL
expr244_mat <- as.matrix(expr244_mat)
mode(expr244_mat) <- "numeric"
meta244 <- read.csv(file.path(root, "analysis/04_adenomyosis/metadata/GSE244236_metadata.csv"), check.names = FALSE)
meta244$sample_id <- sub(".*\\[([^]]+)\\].*", "\\1", meta244$title)
meta244$group <- ifelse(grepl("control", meta244$title, ignore.case = TRUE), "Control", "Adenomyosis")
common244 <- intersect(colnames(expr244_mat), meta244$sample_id)
expr244_mat <- expr244_mat[, common244, drop = FALSE]
meta244 <- meta244[match(common244, meta244$sample_id), ]
common_adeno <- rowMeans(cbind(
  scale_in_dataset(-compute_signature_score(expr244_mat, late_genes, early_genes)),
  scale_in_dataset(compute_signature_score(expr244_mat, rif_up, rif_down)),
  scale_in_dataset(-compute_signature_score(expr244_mat, decidual_genes, character()))
), na.rm = TRUE)
design_adeno <- model.matrix(~ common_adeno + I(meta244$group == "Adenomyosis"))
colnames(design_adeno)[3] <- "case"
fit_adeno <- eBayes(lmFit(expr244_mat, design_adeno))
tt_adeno <- topTable(fit_adeno, coef = "case", number = Inf, sort.by = "none")
stats_adeno <- tt_adeno$t
names(stats_adeno) <- rownames(tt_adeno)
gene_level_residual_work[["Adenomyosis_GO"]] <- run_gsego_table(stats_adeno, "Adenomyosis residual deviation")
if (!is.null(hallmark_gmt)) hallmark_work[["Adenomyosis_H"]] <- run_fgsea_gmt(stats_adeno, hallmark_gmt, "Adenomyosis residual deviation")

gsea_all <- bind_rows(gene_level_residual_work)
write.csv(gsea_all, file.path(tabdir, "disease_specific_residual_fgsea.csv"), row.names = FALSE)

hallmark_all <- bind_rows(hallmark_work)
if (nrow(hallmark_all) > 0) {
  hallmark_all <- hallmark_all %>%
    mutate(across(where(is.list), ~ vapply(.x, function(v) paste(v, collapse = ";"), character(1))))
}
write.csv(hallmark_all, file.path(tabdir, "disease_specific_residual_hallmark.csv"), row.names = FALSE)

gsea_plot_tbl <- gsea_all %>%
  group_by(contrast) %>%
  arrange(p.adjust, desc(abs(NES))) %>%
  mutate(direction = ifelse(NES > 0, "positive", "negative")) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    disease = recode(
      contrast,
      "RIF residual deviation" = "RIF",
      "Endometriosis residual deviation" = "Endometriosis",
      "Adenomyosis residual deviation" = "Adenomyosis"
    ),
    pathway = wrap_text(Description, width = 32)
  )

p_gsea <- ggplot(gsea_plot_tbl, aes(x = NES, y = reorder(pathway, NES), color = NES > 0)) +
  geom_segment(aes(x = 0, xend = NES, yend = reorder(pathway, NES)), linewidth = 0.9, alpha = 0.75) +
  geom_point(aes(size = -log10(p.adjust)), alpha = 0.95) +
  facet_wrap(~ disease, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"), guide = "none") +
  scale_size_continuous(range = c(2.5, 6)) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "#f0f4f8", color = "#d0d7de"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Disease-specific deviation pathways after removing the shared failure axis",
    x = "Residual pathway NES",
    y = "",
    size = "-log10(FDR)",
    subtitle = "Representative Gene Ontology biological-process enrichments from residualized disease contrasts"
  )
ggsave(file.path(figdir, "Figure_ATLAS_RS_10_disease_specific_residual_gsea.png"), p_gsea, width = 10.2, height = 8.2, dpi = 260)

if (nrow(hallmark_all) > 0) {
  hallmark_plot_tbl <- hallmark_all %>%
    group_by(contrast) %>%
    arrange(p.adjust, desc(abs(NES))) %>%
    slice_head(n = 4) %>%
    ungroup() %>%
    mutate(
      disease = recode(
        contrast,
        "RIF residual deviation" = "RIF",
        "Endometriosis residual deviation" = "Endometriosis",
        "Adenomyosis residual deviation" = "Adenomyosis"
      ),
      pathway = wrap_text(gsub("^HALLMARK_", "", Description), width = 26)
    )

  p_hallmark <- ggplot(hallmark_plot_tbl, aes(x = NES, y = reorder(pathway, NES), color = NES > 0)) +
    geom_segment(aes(x = 0, xend = NES, yend = reorder(pathway, NES)), linewidth = 0.9, alpha = 0.75) +
    geom_point(aes(size = -log10(p.adjust)), alpha = 0.95) +
    facet_wrap(~ disease, scales = "free_y", ncol = 1) +
    scale_color_manual(values = c("TRUE" = "#b2182b", "FALSE" = "#2166ac"), guide = "none") +
    scale_size_continuous(range = c(2.5, 6)) +
    theme_bw(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "#f0f4f8", color = "#d0d7de"),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "Hallmark validation of disease-specific residual programs",
      subtitle = "Hallmark pathway enrichments using MSigDB Hallmark gene sets (version 2026.1)",
      x = "Residual pathway NES",
      y = "",
      size = "-log10(FDR)"
    )
  ggsave(file.path(figdir, "Figure_ATLAS_RS_13_hallmark_validation.png"), p_hallmark, width = 10.2, height = 7.8, dpi = 260)
}

# Spatial niche program enrichment network
program_sets <- list(
  CXCL_ligands = c("CXCL1","CXCL2","CXCL3","CXCL8","CXCL9","CXCL10","CXCL11","CXCL12","CXCL13","CXCL14","CXCL16"),
  CCL_ligands = c("CCL2","CCL3","CCL4","CCL5","CCL7","CCL8","CCL19","CCL20","CCL21"),
  IFN_program = c("IFIT1","IFIT2","IFIT3","ISG15","MX1","STAT1","IRF7","CXCL10","CXCL11","GBP1","IFI6"),
  TNF_IL1_program = c("TNF","TNFRSF1A","TNFAIP3","NFKBIA","IL1B","IL1R1","IL1RN","IRAK2","ICAM1"),
  TGFB_VEGF_axis = c("TGFB1","TGFBR1","TGFBR2","SMAD3","VEGFA","FLT1","KDR","ENG"),
  Chemokine_receptors = c("CXCR3","CXCR4","CXCR5","CCR1","CCR2","CCR5","CCR7"),
  MMP_ECM_axis = c("MMP2","MMP7","MMP9","COL1A1","COL3A1","FN1","ITGA5","ITGB1","SPP1")
)

spatial_scores <- read.csv(file.path(root, "analysis/03_rif/round5_maximal/tables/GSE287278_spatial_signature_scores_round5.csv"), check.names = FALSE)
outer_tar <- "/Volumes/Extreme SSD/03_rif/GSE287278/GSE287278_RAW.tar"
outer_members <- utils::untar(outer_tar, list = TRUE)
sample_members <- outer_members[grepl("_processed_data\\.tar\\.gz$", outer_members)]
tmp_outer <- tempfile("gse287278_outer_")
dir.create(tmp_outer)
tmp_inner <- tempfile("gse287278_inner_")
dir.create(tmp_inner)
selected_genes <- unique(unlist(program_sets))
spatial_program_list <- list()

for (member in sample_members) {
  utils::untar(outer_tar, exdir = tmp_outer, files = member)
  inner_path <- file.path(tmp_outer, member)
  sample_stub <- sub("^GSM[0-9]+_", "", basename(member))
  sample_name <- sub("_processed_data\\.tar\\.gz$", "", sample_stub)
  sample_dir <- file.path(tmp_inner, sample_name)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  utils::untar(inner_path, exdir = sample_dir)
  inner_base <- list.files(sample_dir, full.names = TRUE)[1]
  feat_path <- file.path(inner_base, "filtered_feature_bc_matrix", "features.tsv.gz")
  bar_path <- file.path(inner_base, "filtered_feature_bc_matrix", "barcodes.tsv.gz")
  mtx_path <- file.path(inner_base, "filtered_feature_bc_matrix", "matrix.mtx.gz")
  feats <- read.delim(gzfile(feat_path), header = FALSE, stringsAsFactors = FALSE)
  genes <- make.unique(as.character(feats[[2]]))
  keep_idx <- which(genes %in% selected_genes)
  bars <- read.delim(gzfile(bar_path), header = FALSE, stringsAsFactors = FALSE)[[1]]
  mtx <- readMM(mtx_path)
  mtx_sub <- mtx[keep_idx, , drop = FALSE]
  rownames(mtx_sub) <- genes[keep_idx]
  colnames(mtx_sub) <- bars
  libsize <- Matrix::colSums(mtx)
  norm <- log1p(t(t(as.matrix(mtx_sub)) / pmax(libsize, 1)) * 1e4)
  prog_df <- data.frame(
    sample = sample_name,
    barcode = colnames(norm),
    stringsAsFactors = FALSE
  )
  for (nm in names(program_sets)) {
    prog_df[[nm]] <- compute_signature_score_raw(norm, program_sets[[nm]], character())
  }
  spatial_program_list[[sample_name]] <- prog_df
}

spatial_program <- bind_rows(spatial_program_list) %>%
  left_join(spatial_scores, by = c("sample", "barcode"))
write.csv(spatial_program, file.path(tabdir, "rif_spatial_program_scores.csv"), row.names = FALSE)

spatial_program <- spatial_program %>%
  group_by(sample) %>%
  mutate(
    immune_q75 = quantile(rif_immune_score, 0.75, na.rm = TRUE),
    receptivity_q25 = quantile(receptivity_score, 0.25, na.rm = TRUE),
    niche_spot = rif_immune_score >= immune_q75 & receptivity_score <= receptivity_q25
  ) %>%
  ungroup()

program_diff <- bind_rows(lapply(names(program_sets), function(nm) {
  spatial_program %>%
    group_by(sample, group) %>%
    summarise(
      program = nm,
      niche_mean = mean(.data[[nm]][niche_spot], na.rm = TRUE),
      other_mean = mean(.data[[nm]][!niche_spot], na.rm = TRUE),
      delta = niche_mean - other_mean,
      cor_immune = suppressWarnings(cor(.data[[nm]], rif_immune_score, use = "pairwise.complete.obs", method = "spearman")),
      cor_receptivity = suppressWarnings(cor(.data[[nm]], receptivity_score, use = "pairwise.complete.obs", method = "spearman")),
      .groups = "drop"
    )
})) 
write.csv(program_diff, file.path(tabdir, "rif_spatial_niche_program_differences.csv"), row.names = FALSE)

program_summary <- program_diff %>%
  group_by(program) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    lower = bootstrap_mean_ci(delta)["lower"],
    upper = bootstrap_mean_ci(delta)["upper"],
    mean_cor_immune = mean(cor_immune, na.rm = TRUE),
    mean_cor_receptivity = mean(cor_receptivity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_delta))
write.csv(program_summary, file.path(tabdir, "rif_spatial_niche_program_summary.csv"), row.names = FALSE)

program_quant <- program_diff %>%
  group_by(program) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    lower = bootstrap_mean_ci(delta)["lower"],
    upper = bootstrap_mean_ci(delta)["upper"],
    n_sections = sum(is.finite(delta)),
    same_direction_sections = sum(sign(delta[is.finite(delta)]) == sign(mean(delta, na.rm = TRUE))),
    signflip_p = exact_signflip_p(delta),
    .groups = "drop"
  ) %>%
  mutate(
    fdr = p.adjust(signflip_p, method = "BH"),
    consistency = paste0(same_direction_sections, "/", n_sections, " sections")
  ) %>%
  arrange(fdr, desc(abs(mean_delta)))
write.csv(program_quant, file.path(tabdir, "rif_spatial_niche_program_quantitative_summary.csv"), row.names = FALSE)

network_df <- program_summary %>%
  mutate(
    family = case_when(
      grepl("CXCL", program) ~ "Chemokine ligand",
      grepl("CCL", program) ~ "Chemokine ligand",
      grepl("receptors", program) ~ "Receptor axis",
      grepl("ECM", program) ~ "ECM / remodeling",
      TRUE ~ "Inflammatory / growth axis"
    ),
    angle = seq(pi / 2, pi / 2 + 2 * pi, length.out = n() + 1)[1:n()],
    x = cos(angle),
    y = sin(angle),
    xend = 0,
    yend = 0
  ) %>%
  arrange(desc(abs(mean_delta))) %>%
  mutate(program_label = wrap_text(gsub("_", " ", program), width = 16))

p_network <- ggplot() +
  geom_segment(
    data = network_df,
    aes(x = x, y = y, xend = xend, yend = yend, linewidth = abs(mean_delta), color = family),
    alpha = 0.75
  ) +
  geom_point(data = network_df, aes(x = x, y = y, fill = mean_cor_immune), shape = 21, color = "black", size = 8.5, stroke = 0.5) +
  geom_text(data = network_df, aes(x = x * 1.34, y = y * 1.34, label = program_label), size = 3.5, lineheight = 0.95) +
  geom_point(aes(x = 0, y = 0), size = 21, shape = 21, fill = "#6a3d9a", color = "black", stroke = 0.7) +
  annotate("text", x = 0, y = 0, label = "Joint\nniche", color = "white", size = 3.1, fontface = "bold", lineheight = 0.92) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", name = "Mean spearman\ncorr. with immune") +
  scale_linewidth_continuous(name = "|Niche delta|", range = c(0.8, 3.2)) +
  coord_equal(xlim = c(-1.82, 1.82), ylim = c(-1.65, 1.72), clip = "off") +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 10.5),
    plot.background = element_rect(fill = "white", color = "white"),
    panel.background = element_rect(fill = "white", color = "white"),
    legend.background = element_rect(fill = "white", color = "white"),
    legend.key = element_rect(fill = "white", color = "white")
  ) +
  labs(
    title = "Spot-level ligand/receptor program enrichment around immune-high and receptivity-low niches",
    subtitle = "Program nodes summarize section-level niche enrichment rather than direct cell-cell inference"
  )
ggsave(file.path(figdir, "Figure_ATLAS_RS_11_spatial_niche_program_network.png"), p_network, width = 9.2, height = 6.8, dpi = 260, bg = "white")

# Clinical WOI explainability figure
canon <- read.csv(file.path(tabdir, "canonical_receptivity_overlap.csv"), check.names = FALSE)
canon_long <- canon %>%
  pivot_longer(cols = c(in_dynamic_genes, in_80_gene_signature), names_to = "set_name", values_to = "present") %>%
  mutate(
    set_name = recode(set_name,
                      in_dynamic_genes = "3,499 dynamic genes",
                      in_80_gene_signature = "80-gene transfer signature")
  )

flow_tbl <- data.frame(
  stage = factor(c("Canonical clinical\nreceptivity markers", "Present in dynamic\ngene framework", "Retained in compact\ntransfer signature"),
                 levels = c("Canonical clinical\nreceptivity markers", "Present in dynamic\ngene framework", "Retained in compact\ntransfer signature")),
  n = c(12, sum(canon$in_dynamic_genes), sum(canon$in_80_gene_signature))
)

p_flow <- ggplot(flow_tbl, aes(stage, n, fill = stage)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = n), vjust = -0.4, size = 4.2) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(size = 9)) +
  labs(title = "Clinical WOI marker coverage versus the broader timing framework", x = "", y = "Number of markers") +
  guides(fill = "none")

p_tile <- ggplot(canon_long, aes(set_name, gene, fill = present)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("TRUE" = "#b2182b", "FALSE" = "grey90")) +
  theme_bw(base_size = 10) +
  labs(title = "Canonical receptivity-marker representation", x = "", y = "", fill = "Present") +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

png(file.path(figdir, "Figure_ATLAS_RS_12_clinical_marker_explainability.png"), width = 2200, height = 1400, res = 240)
grid::grid.newpage()
pushViewport <- grid::pushViewport
viewport <- grid::viewport
grid::pushViewport(viewport(layout = grid::grid.layout(1, 2)))
print(p_flow, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_tile, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()