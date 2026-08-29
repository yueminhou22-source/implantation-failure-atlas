#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(lmerTest)
  library(lme4)
})

set.seed(123)

root <- "path omitted"
outdir <- file.path(root, "analysis/05_implantation_failure_atlas/reviewer_stats")
figdir <- file.path(outdir, "figures")
tabdir <- file.path(outdir, "tables")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)

heca_dir <- file.path(root, "analysis/05_implantation_failure_atlas/heca_upgrade/tables")
atlas_path <- file.path(root, "analysis/05_implantation_failure_atlas/final/tables/cross_disease_score_table.csv")
# The manuscript's 3,499-gene definition is the round5 Spearman/BH list,
# not the broader stage-F-statistic table.
dyn_path <- file.path(root, "analysis/03_rif/round5_maximal/tables/fertile_day_spearman_round5.csv")
sig_path <- file.path(root, "analysis/03_rif/round5_maximal/tables/receptivity_signature_genes_round5.csv")

build_sample_metrics <- function(cache_file, dataset_name, phase_map = NULL) {
  obj <- readRDS(cache_file)
  cells <- obj$cells %>%
    group_by(sample, group) %>%
    summarise(
      total_cells = sum(cell_n, na.rm = TRUE),
      n_celltypes = n_distinct(assigned_celltype),
      .groups = "drop"
    )
  agg <- obj$agg %>%
    group_by(sample, group) %>%
    summarise(
      total_counts = sum(value, na.rm = TRUE),
      mean_n_cells_state = mean(n_cells, na.rm = TRUE),
      .groups = "drop"
    )
  summ <- obj$summary %>%
    mutate(dataset = dataset_name)
  out <- summ %>%
    left_join(cells, by = c("sample", "group")) %>%
    left_join(agg, by = c("sample", "group"))
  if (!is.null(phase_map)) {
    out <- out %>% left_join(phase_map, by = "sample")
  }
  out
}

phase_map <- read.csv(file.path(heca_dir, "heca_assignment_summary.csv"), check.names = FALSE) %>%
  select(sample, phase) %>%
  distinct()

sample_metrics <- bind_rows(
  build_sample_metrics(file.path(heca_dir, "rif_heca_raw_cache.rds"), "GSE250130", phase_map),
  build_sample_metrics(file.path(heca_dir, "endo179640_heca_raw_cache.rds"), "GSE179640"),
  build_sample_metrics(file.path(heca_dir, "endo213216_heca_raw_cache.rds"), "GSE213216"),
  build_sample_metrics(file.path(heca_dir, "endo214411_heca_raw_cache.rds"), "GSE214411", phase_map)
) %>%
  mutate(
    dataset = factor(dataset, levels = c("GSE250130", "GSE179640", "GSE213216", "GSE214411")),
    context_class = case_when(
      dataset == "GSE250130" & group == "RIF" ~ "RIF_case",
      dataset == "GSE250130" ~ "fertile_reference",
      dataset == "GSE179640" & group == "Control_Eutopic" ~ "control",
      dataset == "GSE179640" & group == "Endo_Eutopic" ~ "eutopic_endometriosis",
      dataset == "GSE179640" ~ "lesion_associated",
      dataset == "GSE213216" & group == "No_endo_detected" ~ "control",
      dataset == "GSE213216" & group == "Eutopic" ~ "eutopic_endometriosis",
      dataset == "GSE213216" ~ "lesion_associated",
      dataset == "GSE214411" & grepl("Control", group, ignore.case = TRUE) ~ "control",
      dataset == "GSE214411" ~ "eutopic_endometriosis",
      TRUE ~ "other"
    ),
    log10_total_cells = log10(total_cells + 1),
    log10_total_counts = log10(total_counts + 1)
  )
write.csv(sample_metrics, file.path(tabdir, "heca_assignment_sample_metrics.csv"), row.names = FALSE)

lm_dataset <- lm(mean_cor ~ dataset, data = sample_metrics)
lm_size <- lm(mean_cor ~ dataset + log10_total_cells + log10_total_counts, data = sample_metrics)
lm_full <- lm(mean_cor ~ dataset + log10_total_cells + log10_total_counts + context_class, data = sample_metrics)

extract_coef <- function(mod, label) {
  sm <- summary(mod)$coefficients
  data.frame(
    model = label,
    term = rownames(sm),
    estimate = sm[, 1],
    std_error = sm[, 2],
    statistic = sm[, 3],
    p_value = sm[, 4],
    row.names = NULL
  )
}

coef_tbl <- bind_rows(
  extract_coef(lm_dataset, "dataset_only"),
  extract_coef(lm_size, "dataset_plus_size_depth"),
  extract_coef(lm_full, "dataset_plus_size_depth_context")
)
write.csv(coef_tbl, file.path(tabdir, "heca_assignment_covariate_models.csv"), row.names = FALSE)

model_comp <- anova(lm_dataset, lm_size, lm_full)
write.csv(
  data.frame(model = rownames(model_comp), model_comp, row.names = NULL),
  file.path(tabdir, "heca_assignment_model_comparison.csv"),
  row.names = FALSE
)

dataset_summary <- sample_metrics %>%
  group_by(dataset) %>%
  summarise(
    n_samples = n(),
    mean_cor = mean(mean_cor, na.rm = TRUE),
    sd_cor = sd(mean_cor, na.rm = TRUE),
    median_cells = median(total_cells, na.rm = TRUE),
    median_counts = median(total_counts, na.rm = TRUE),
    prop_secretory = mean(phase == "secretory", na.rm = TRUE),
    .groups = "drop"
  )
write.csv(dataset_summary, file.path(tabdir, "heca_assignment_dataset_summary.csv"), row.names = FALSE)

plot_df <- sample_metrics %>%
  select(sample, dataset, mean_cor, context_class, log10_total_cells, log10_total_counts) %>%
  pivot_longer(
    cols = c(log10_total_cells, log10_total_counts),
    names_to = "metric",
    values_to = "metric_value"
  ) %>%
  mutate(
    metric = recode(metric,
                    log10_total_cells = "log10(total cells)",
                    log10_total_counts = "log10(total counts)")
  )

p_cov <- ggplot(plot_df, aes(metric_value, mean_cor, color = dataset, shape = context_class)) +
  geom_point(size = 2.5, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, alpha = 0.3) +
  facet_wrap(~ metric, scales = "free_x") +
  theme_bw(base_size = 12) +
  labs(
    title = "HECA assignment-correlation sensitivity to sample size and sequencing depth",
    x = "",
    y = "Mean HECA assignment correlation",
    color = "Dataset",
    shape = "Context"
  )
ggsave(file.path(figdir, "Figure_ATLAS_RS_6_heca_covariate_sensitivity.png"), p_cov, width = 9.4, height = 5.2, dpi = 260)

calc_shift <- function(df, case_group, ctrl_group) {
  case <- df$pred_day[df$group == case_group]
  ctrl <- df$pred_day[df$group == ctrl_group]
  data.frame(
    n_case = sum(df$group == case_group),
    n_ctrl = sum(df$group == ctrl_group),
    mean_case = mean(case, na.rm = TRUE),
    mean_ctrl = mean(ctrl, na.rm = TRUE),
    delta = mean(case, na.rm = TRUE) - mean(ctrl, na.rm = TRUE)
  )
}

threshold_grid <- c(0, 0.45, 0.50, 0.55)
assign_tbl <- read.csv(file.path(heca_dir, "heca_assignment_summary.csv"), check.names = FALSE)

timing_files <- list(
  rif = list(
    file = file.path(heca_dir, "rif_heca_celltype_timing_scores.csv"),
    dataset = "GSE250130",
    celltype = "Glandular_Epi",
    case_group = "RIF",
    ctrl_group = "Fertile_LH7",
    label = "RIF glandular epithelium"
  ),
  endo179640 = list(
    file = file.path(heca_dir, "endo179640_heca_celltype_timing_scores.csv"),
    dataset = "GSE179640",
    celltype = "Stroma",
    case_group = "Endo_Eutopic",
    ctrl_group = "Control_Eutopic",
    label = "GSE179640 eutopic stroma"
  ),
  endo213216 = list(
    file = file.path(heca_dir, "endo213216_heca_celltype_timing_scores.csv"),
    dataset = "GSE213216",
    celltype = "Decidual_Stroma",
    case_group = "Lesion",
    ctrl_group = "Eutopic",
    label = "GSE213216 lesion decidual stroma"
  ),
  endo214411 = list(
    file = file.path(heca_dir, "endo214411_heca_celltype_timing_scores.csv"),
    dataset = "GSE214411",
    celltype = "Endothelial",
    case_group = "endometriosis",
    ctrl_group = "Control",
    label = "GSE214411 endothelial"
  )
)

threshold_sensitivity <- bind_rows(lapply(timing_files, function(spec) {
  df <- read.csv(spec$file, check.names = FALSE) %>%
    filter(celltype == spec$celltype, group %in% c(spec$case_group, spec$ctrl_group)) %>%
    left_join(
      assign_tbl %>% filter(dataset == spec$dataset) %>% select(sample, mean_cor),
      by = "sample"
    )
  bind_rows(lapply(threshold_grid, function(thr) {
    out <- calc_shift(df %>% filter(mean_cor >= thr), spec$case_group, spec$ctrl_group)
    cbind(
      contrast = spec$label,
      dataset = spec$dataset,
      celltype = spec$celltype,
      threshold = thr,
      out
    )
  }))
}))
write.csv(threshold_sensitivity, file.path(tabdir, "heca_assignment_threshold_sensitivity.csv"), row.names = FALSE)

atlas <- read.csv(atlas_path, check.names = FALSE)

case_df <- atlas %>%
  filter(!(dataset %in% c("GSE287278", "GSE263897"))) %>%
  filter(
    (disease == "RIF" & group == "RIF") |
      (disease == "Endometriosis" & !(group %in% c("Control_Eutopic", "Control", "healthy control", "No_endo_detected"))) |
      (disease == "Adenomyosis" & group == "Adenomyosis")
  ) %>%
  mutate(
    disease = factor(disease, levels = c("RIF", "Endometriosis", "Adenomyosis")),
    material_class = case_when(
      disease == "RIF" ~ "endometrium_like",
      disease == "Endometriosis" & group %in% c("Endo_Eutopic", "Eutopic", "endometriosis", "patient with endometriosis") ~ "endometrium_like",
      disease == "Endometriosis" ~ "lesion_like",
      disease == "Adenomyosis" & dataset == "GSE244236" ~ "organoid",
      TRUE ~ "tissue_or_stroma"
    )
  )

case_model <- lmer(pred_day ~ disease + (1 | dataset), data = case_df, REML = FALSE)
case_coef <- coef(summary(case_model))
case_coef_out <- data.frame(
  term = rownames(case_coef),
  estimate = case_coef[, "Estimate"],
  std_error = case_coef[, "Std. Error"],
  df = case_coef[, "df"],
  statistic = case_coef[, "t value"],
  p_value = case_coef[, "Pr(>|t|)"],
  row.names = NULL
)
write.csv(case_coef_out, file.path(tabdir, "disease_case_mixed_model_coefficients.csv"), row.names = FALSE)

case_anova <- anova(case_model)
write.csv(
  data.frame(term = rownames(case_anova), case_anova, row.names = NULL),
  file.path(tabdir, "disease_case_mixed_model_anova.csv"),
  row.names = FALSE
)

pairwise_from_ref <- function(ref_level) {
  dat <- case_df
  dat$disease <- relevel(dat$disease, ref = ref_level)
  mod <- lmer(pred_day ~ disease + (1 | dataset), data = dat, REML = FALSE)
  cf <- coef(summary(mod))
  out <- data.frame(
    reference = ref_level,
    contrast = rownames(cf),
    estimate = cf[, "Estimate"],
    std_error = cf[, "Std. Error"],
    df = cf[, "df"],
    statistic = cf[, "t value"],
    p_value = cf[, "Pr(>|t|)"],
    row.names = NULL
  ) %>%
    filter(contrast != "(Intercept)")
  out
}

pairwise_tbl <- bind_rows(
  pairwise_from_ref("RIF"),
  pairwise_from_ref("Endometriosis"),
  pairwise_from_ref("Adenomyosis")
) %>%
  mutate(
    lower = estimate - 1.96 * std_error,
    upper = estimate + 1.96 * std_error
  )
write.csv(pairwise_tbl, file.path(tabdir, "disease_case_pairwise_contrasts.csv"), row.names = FALSE)

disease_summary <- case_df %>%
  group_by(disease) %>%
  summarise(
    n = n(),
    mean_pred_day = mean(pred_day, na.rm = TRUE),
    lower = quantile(replicate(2000, mean(sample(pred_day, replace = TRUE), na.rm = TRUE)), 0.025, na.rm = TRUE),
    upper = quantile(replicate(2000, mean(sample(pred_day, replace = TRUE), na.rm = TRUE)), 0.975, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(disease_summary, file.path(tabdir, "disease_case_pred_day_summary.csv"), row.names = FALSE)

p_disease <- ggplot(case_df, aes(disease, pred_day, color = disease)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1.5) +
  geom_pointrange(
    data = disease_summary,
    aes(x = disease, y = mean_pred_day, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.7
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Predicted receptive day differs across disease-case observations",
    subtitle = "Black points and bars indicate bootstrap mean and 95% CI; model treats dataset as a random intercept",
    x = "",
    y = "Predicted receptive day"
  ) +
  guides(color = "none")
ggsave(file.path(figdir, "Figure_ATLAS_RS_7_disease_case_pred_day.png"), p_disease, width = 7.6, height = 5.2, dpi = 260)

dyn <- read.csv(dyn_path, check.names = FALSE)
sig <- read.csv(sig_path, check.names = FALSE)
dyn$gene <- toupper(trimws(dyn$gene))
dynamic_genes <- dyn$gene[dyn$padj < 0.05]
canonical <- data.frame(
  gene = c("PAEP", "GPX3", "MUC1", "IGFBP1", "PRL", "LEFTY2", "SPP1", "HOXA10", "HOXA11", "IL15", "CXCL14", "DPP4"),
  category = "canonical_receptivity_marker"
)
sig_genes <- unique(c(na.omit(sig$late_genes), na.omit(sig$early_genes)))
canonical_overlap <- canonical %>%
  mutate(
    in_dynamic_genes = gene %in% dynamic_genes,
    in_80_gene_signature = gene %in% sig_genes
  )
write.csv(canonical_overlap, file.path(tabdir, "canonical_receptivity_overlap.csv"), row.names = FALSE)

summary_overlap <- data.frame(
  gene_set = c("Canonical ERA-like clinical receptivity markers"),
  n_total = nrow(canonical_overlap),
  n_in_dynamic_genes = sum(canonical_overlap$in_dynamic_genes),
  pct_in_dynamic_genes = 100 * mean(canonical_overlap$in_dynamic_genes),
  n_in_80_gene_signature = sum(canonical_overlap$in_80_gene_signature),
  pct_in_80_gene_signature = 100 * mean(canonical_overlap$in_80_gene_signature)
)
write.csv(summary_overlap, file.path(tabdir, "canonical_receptivity_overlap_summary.csv"), row.names = FALSE)