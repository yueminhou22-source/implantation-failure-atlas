suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

root <- getwd()
source_dir <- file.path(root, "source_data")
output_dir <- file.path(root, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

theme_submission <- theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 18, 10, 10, "pt")
  )

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(output_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 600, limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(stem, ".tiff")), plot, width = width, height = height, dpi = 600, compression = "lzw", limitsize = FALSE)
  ggsave(file.path(output_dir, paste0(stem, ".pdf")), plot, width = width, height = height, limitsize = FALSE)
}

# Figure S1: identical source data and geometry; only the display legend title and labels are publication-facing.
loocv <- read_csv(file.path(source_dir, "fertile_timeline_loocv_predictions.csv"), show_col_types = FALSE)
metrics <- read_csv(file.path(source_dir, "fertile_timeline_loocv_metrics.csv"), show_col_types = FALSE)
metric_value <- function(name) metrics$value[metrics$metric == name][1]
loocv$stage_display <- factor(
  loocv$group_simple,
  levels = c("LH3", "LH5", "Fertile_LH7", "LH9", "LH11"),
  labels = c("LH3", "LH5", "Fertile LH7", "LH9", "LH11")
)
p_s1 <- ggplot(loocv, aes(true_day, pred_day, color = stage_display)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 3, alpha = 0.9, position = position_jitter(width = 0.08, height = 0.08)) +
  scale_x_continuous(breaks = c(3, 5, 7, 9, 11)) +
  scale_y_continuous(breaks = c(3, 5, 7, 9, 11)) +
  scale_color_manual(
    name = "Luteal stage",
    values = c("LH3" = "#2166ac", "LH5" = "#67a9cf", "Fertile LH7" = "#1b7837", "LH9" = "#ef8a62", "LH11" = "#b2182b")
  ) +
  annotate(
    "text", x = 10.8, y = 4.15, hjust = 1,
    label = paste0(
      "LOOCV rho = ", sprintf("%.2f", metric_value("LOOCV Spearman rho")), "\n",
      "MAE = ", sprintf("%.2f", metric_value("LOOCV mean absolute error")), " days\n",
      "Exact stage = ", sprintf("%.0f%%", 100 * metric_value("LOOCV exact-stage recovery")), "\n",
      "Pairwise order = ", sprintf("%.0f%%", 100 * metric_value("LOOCV pairwise stage-order concordance"))
    ), size = 3.8, lineheight = 1.05
  ) +
  theme_submission +
  labs(title = "Leave-one-sample-out validation of the fertile receptivity timeline", x = "Observed luteal day", y = "Predicted luteal day")
save_plot(p_s1, "Figure_S1_rebuilt_display_labels", 7.4, 5.8)

# Figure S12: source categories are author-annotated/fine categories, not a substituted HECA broad-state ontology.
s12 <- read_csv(file.path(source_dir, "rif_celltype_timing_stats.csv"), show_col_types = FALSE) %>%
  mutate(
    display_celltype = recode(
      celltype,
      "Decidual_Stroma" = "Decidual stroma",
      "Endothelial" = "Endothelium",
      "Glandular_Epi" = "Glandular epithelium",
      "Luminal_Epi" = "Luminal epithelium",
      "Macrophage" = "Macrophages",
      "NK_T" = "NK/T cells",
      "Stroma" = "Stroma"
    ),
    display_celltype = factor(display_celltype, levels = display_celltype[order(estimate)])
  )
p_s12 <- ggplot(s12, aes(estimate, display_celltype, color = estimate > 0)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18, linewidth = 0.7) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c("FALSE" = "#2C7FB8", "TRUE" = "#D95F0E"), guide = "none") +
  theme_submission +
  labs(title = "Cell-type-resolved timing displacement in RIF", x = "Predicted-day shift in RIF versus fertile LH7", y = "")
save_plot(p_s12, "Figure_S12_rebuilt_display_labels", 7.6, 5.4)

# Figure S14: source counts are HECA reference coverage, not disease effects.
s14 <- read_csv(file.path(source_dir, "heca_broad_celltype_stage_counts.csv"), show_col_types = FALSE) %>%
  filter(stage_simple %in% c("Secretory_Early", "Secretory_EarlyMid", "Secretory_Mid", "Secretory_Late")) %>%
  mutate(
    stage_display = factor(
      stage_simple,
      levels = c("Secretory_Early", "Secretory_EarlyMid", "Secretory_Mid", "Secretory_Late"),
      labels = c("Secretory early", "Secretory early-mid", "Secretory mid", "Secretory late")
    ),
    celltype_display = recode(
      broad_celltype,
      "Decidual_Stroma" = "Decidual stroma",
      "Endothelial" = "Endothelium",
      "Glandular_Epi" = "Glandular epithelium",
      "Luminal_Epi" = "Luminal epithelium",
      "Stroma" = "Stroma",
      "Lymphoid" = "Lymphoid",
      "Myeloid" = "Myeloid"
    )
  )
p_s14 <- ggplot(s14, aes(stage_display, celltype_display, fill = n_cells)) +
  geom_tile() +
  scale_fill_gradient(low = "#F4E8C1", high = "#A73A24", name = "HECA reference\ncells") +
  theme_submission +
  labs(title = "HECA secretory-stage coverage across broad endometrial cell types", x = "", y = "")
save_plot(p_s14, "Figure_S14_rebuilt_display_labels", 8.2, 4.8)

# Figure S15: audit and rebuild at comparable sample/pseudobulk/organoid analytical units only.
s15 <- read_csv(file.path(source_dir, "cross_disease_score_table_with_common_axis.csv"), show_col_types = FALSE)
unit_map <- c(
  "GSE111974" = "bulk sample", "GSE58144" = "bulk sample", "GSE135485" = "bulk sample",
  "GSE157718" = "bulk sample", "GSE190580" = "bulk sample", "GSE78851" = "bulk sample",
  "GSE179640" = "scRNA pseudobulk sample", "GSE213216" = "scRNA pseudobulk sample",
  "GSE214411" = "scRNA pseudobulk sample", "GSE250130" = "scRNA pseudobulk sample",
  "GSE287278" = "Visium spot", "GSE263897" = "GeoMx ROI", "GSE244236" = "organoid sample"
)
type_map <- c(
  "GSE111974" = "Bulk transcriptomics", "GSE58144" = "Bulk transcriptomics", "GSE135485" = "Bulk transcriptomics",
  "GSE157718" = "Stromal transcriptomics", "GSE190580" = "Tissue transcriptomics", "GSE78851" = "Whole-tissue transcriptomics",
  "GSE179640" = "Single-cell-derived pseudobulk", "GSE213216" = "Single-cell-derived pseudobulk",
  "GSE214411" = "Single-cell-derived pseudobulk", "GSE250130" = "Single-cell-derived pseudobulk",
  "GSE287278" = "Visium spatial transcriptomics", "GSE263897" = "GeoMx spatial transcriptomics", "GSE244236" = "Organoid transcriptomics"
)
case_flag <- with(s15,
  (disease == "RIF" & group == "RIF") |
  (disease == "Endometriosis" & !(group %in% c("Control_Eutopic", "Control", "healthy control", "No_endo_detected"))) |
  (disease == "Adenomyosis" & group == "Adenomyosis")
)
s15_audit <- s15 %>%
  mutate(
    observation_id = sample,
    original_data_type = unname(type_map[dataset]),
    original_inference_unit = unname(unit_map[dataset]),
    plot_unit = original_inference_unit,
    is_case_final lock = case_flag,
    included_in_legacy_S15 = is_case_final lock & !is.na(implantation_module_score),
    included_in_S15 = included_in_legacy_S15 & original_inference_unit %in% c("bulk sample", "scRNA pseudobulk sample", "organoid sample"),
    reason_for_inclusion = case_when(
      !is_case_final lock ~ "Not a disease-case observation under the original Figure S15 case definition",
      is.na(implantation_module_score) ~ "Excluded because fewer than two implantation-module components were available",
      original_inference_unit == "Visium spot" ~ "Excluded from cleaned Figure S15 because Visium spots are spatial analytical observations",
      original_inference_unit == "GeoMx ROI" ~ "Excluded from cleaned Figure S15 because GeoMx ROIs are spatial analytical observations",
      TRUE ~ "Retained as a comparable bulk, sample pseudobulk, or organoid analytical observation"
    )
  ) %>%
  group_by(dataset) %>%
  mutate(n_observations = n()) %>%
  ungroup() %>%
  select(dataset, disease, observation_id, original_data_type, original_inference_unit, plot_unit, n_observations, included_in_legacy_S15, included_in_S15, reason_for_inclusion, implantation_module_score)
write_csv(s15_audit, file.path(source_dir, "figure_S15_observation_unit_audit.csv"))

s15_summary <- s15_audit %>%
  filter(included_in_legacy_S15) %>%
  count(dataset, disease, original_data_type, original_inference_unit, plot_unit, included_in_S15, name = "n_observations")
write_csv(s15_summary, file.path(source_dir, "figure_S15_observation_unit_summary.csv"))

s15_clean <- s15_audit %>%
  filter(included_in_S15) %>%
  mutate(disease = factor(disease, levels = c("Adenomyosis", "Endometriosis", "RIF")))
p_s15 <- ggplot(s15_clean, aes(disease, implantation_module_score, color = disease)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.8) +
  geom_jitter(width = 0.14, alpha = 0.48, size = 1.55) +
  scale_color_manual(values = c("RIF" = "#619CFF", "Endometriosis" = "#00BA38", "Adenomyosis" = "#F8766D"), guide = "none") +
  theme_submission +
  labs(
    title = "Hypothesis-driven implantation-failure module score across disorders",
    x = "",
    y = "Implantation-failure module score"
  )
save_plot(p_s15, "Figure_S15_rebuilt_sample_level", 7.2, 5.2)

cat("Rebuilt S1, S12, S14, and S15 with audited source data and publication-facing labels.\n")