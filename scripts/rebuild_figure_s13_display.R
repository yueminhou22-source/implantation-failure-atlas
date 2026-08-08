#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

root <- normalizePath(".", mustWork = TRUE)
source_dir <- file.path(root, "source_data")
out_dir <- file.path(root, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

display_celltype <- function(x) {
  dplyr::recode(
    x,
    "Decidual_Stroma" = "Decidual stroma",
    "Glandular_Epi" = "Glandular epithelium",
    "Luminal_Epi" = "Luminal epithelium",
    "Endothelial" = "Endothelial",
    "Lymphoid" = "Lymphoid",
    "Myeloid" = "Myeloid",
    "Stroma" = "Stroma",
    .default = x
  )
}

display_group <- function(x) {
  dplyr::recode(
    x,
    "Control_Eutopic" = "Control eutopic",
    "Endo_Eutopic" = "Endometriosis eutopic",
    "Ectopic_Ovarian" = "Ectopic ovarian lesion",
    "Ectopic_Adjacent" = "Ectopic adjacent",
    "No_endo_detected" = "No endometriosis detected",
    "Fertile_LH7" = "Fertile LH7",
    "RIF" = "RIF",
    "Control" = "Control",
    "endometriosis" = "Endometriosis",
    "Lesion" = "Lesion",
    .default = x
  )
}

load_scores <- function(file, dataset_label, rif = FALSE) {
  x <- read_csv(file.path(source_dir, file), show_col_types = FALSE)
  if (rif) {
    x <- x %>% mutate(group = case_when(
      group == "RIF" ~ "RIF",
      grepl("LH7", sample) ~ "Fertile_LH7",
      TRUE ~ as.character(group)
    ))
  }
  x %>%
    mutate(dataset = dataset_label,
           celltype_display = display_celltype(celltype),
           group_display = display_group(group),
           context = paste(dataset, group_display, sep = " | ")) %>%
    group_by(context, celltype_display) %>%
    summarise(mean_pred_day = mean(pred_day, na.rm = TRUE), .groups = "drop")
}

heat_data <- bind_rows(
  load_scores("rif_heca_celltype_timing_scores.csv", "GSE250130", rif = TRUE),
  load_scores("endo179640_heca_celltype_timing_scores.csv", "GSE179640"),
  load_scores("endo213216_heca_celltype_timing_scores.csv", "GSE213216"),
  load_scores("endo214411_heca_celltype_timing_scores.csv", "GSE214411")
) %>%
  filter(celltype_display != "Unknown")

# Preserve first-appearance context ordering and a biological broad-state ordering.
row_order <- c("Luminal epithelium", "Glandular epithelium", "Stroma", "Decidual stroma", "Endothelial", "Lymphoid", "Myeloid")
col_order <- unique(heat_data$context)
plot_data <- heat_data %>%
  mutate(
    celltype_display = factor(celltype_display, levels = rev(row_order)),
    context = factor(context, levels = col_order)
  )

write_csv(heat_data, file.path(source_dir, "Figure_S13_cross_cohort_timing_source.csv"))

plot <- ggplot(plot_data, aes(context, celltype_display, fill = mean_pred_day)) +
  geom_tile(color = "white", linewidth = 0.45) +
  geom_text(aes(label = sprintf("%.1f", mean_pred_day)), size = 3.0, color = "black") +
  scale_fill_gradient2(low = "#3B4CC0", mid = "#F7F7F7", high = "#B40426", midpoint = 7, name = "Mean predicted\nreceptive day") +
  labs(x = NULL, y = NULL) +
  # Helvetica is available to the local R graphics devices and is visually
  # consistent with the journal submission font specification.
  theme_minimal(base_size = 11, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(12, 22, 28, 12)
  )

ggsave(file.path(out_dir, "Figure_S13_rebuilt_display_labels.png"), plot, width = 10.2, height = 5.5, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Figure_S13_rebuilt_display_labels.tiff"), plot, width = 10.2, height = 5.5, dpi = 600, compression = "lzw", limitsize = FALSE)
ggsave(file.path(out_dir, "Figure_S13_rebuilt_display_labels.pdf"), plot, width = 10.2, height = 5.5, limitsize = FALSE)