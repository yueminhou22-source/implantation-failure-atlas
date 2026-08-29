#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: rebuild_corrected_figure1c_s1_s3.R <corrected_root> <output_root>")
}
corrected_root <- normalizePath(args[[1]], mustWork = TRUE)
output_root <- normalizePath(args[[2]], mustWork = TRUE)

dir.create(file.path(output_root, "01_Figure1C"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "02_Supplementary_S1"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "03_Supplementary_S3"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_root, "05_logs"), recursive = TRUE, showWarnings = FALSE)

primary_dir <- file.path(corrected_root, "02_corrected_primary_loocv")
null_dir <- file.path(corrected_root, "03_corrected_negative_controls")
pred_path <- file.path(primary_dir, "corrected_primary_loocv_predictions.csv")
summary_path <- file.path(primary_dir, "corrected_primary_loocv_summary.csv")
fold_audit_path <- file.path(primary_dir, "corrected_primary_loocv_fold_audit.csv")
random_path <- file.path(null_dir, "corrected_random80_null_metrics.csv")
permuted_path <- file.path(null_dir, "corrected_permutedLH_null_metrics.csv")
low_assoc_path <- file.path(null_dir, "corrected_low_association_null_metrics.csv")
null_summary_path <- file.path(null_dir, "corrected_negative_control_summary.csv")

required <- c(pred_path, summary_path, fold_audit_path, random_path, permuted_path,
              low_assoc_path, null_summary_path)
if (!all(file.exists(required))) {
  stop("Missing corrected input(s): ", paste(required[!file.exists(required)], collapse = "; "))
}

primary <- read.csv(pred_path, check.names = FALSE, stringsAsFactors = FALSE)
summary <- read.csv(summary_path, check.names = FALSE, stringsAsFactors = FALSE)
random_null <- read.csv(random_path, check.names = FALSE, stringsAsFactors = FALSE)
permuted_null <- read.csv(permuted_path, check.names = FALSE, stringsAsFactors = FALSE)
low_assoc_null <- read.csv(low_assoc_path, check.names = FALSE, stringsAsFactors = FALSE)
null_summary <- read.csv(null_summary_path, check.names = FALSE, stringsAsFactors = FALSE)

if (nrow(primary) != 18) stop("Corrected primary prediction count is not 18")
if (nrow(random_null) != 250 || nrow(permuted_null) != 250 || nrow(low_assoc_null) != 250) {
  stop("Corrected null-control iteration count is not 250 for all controls")
}

rho <- as.numeric(summary$rho[[1]])
mae <- as.numeric(summary$mae[[1]])
pairwise <- as.numeric(summary$pairwise_concordance[[1]])
exact <- as.numeric(summary$exact_nearest_stage_agreement[[1]])

stage_labels <- c(`3` = "LH3", `5` = "LH5", `7` = "Fertile_LH7", `9` = "LH9", `11` = "LH11")
primary$group_simple <- unname(stage_labels[as.character(primary$observed_day)])
primary$group_simple <- factor(primary$group_simple,
                               levels = c("Fertile_LH7", "LH11", "LH3", "LH5", "LH9"))

palette <- c(Fertile_LH7 = "#2E8B57", LH11 = "#C8233B", LH3 = "#2F75B5",
             LH5 = "#6BAED6", LH9 = "#F28E68")

theme_figure1c <- theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(
    panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.45),
    panel.grid.minor = element_line(colour = "#F0F0F0", linewidth = 0.35),
    panel.border = element_rect(colour = "#555555", linewidth = 0.6, fill = NA),
    axis.line = element_line(colour = "#303030", linewidth = 0.45),
    axis.ticks = element_line(colour = "#303030", linewidth = 0.45),
    axis.text = element_text(size = 10, colour = "#303030"),
    axis.title = element_text(size = 12, colour = "#202020"),
    legend.title = element_text(size = 10.5, colour = "#303030"),
    legend.text = element_text(size = 10, colour = "#303030"),
    legend.key = element_blank(),
    plot.margin = margin(8, 18, 8, 10, unit = "pt")
  )

set.seed(20260824)
figure1c <- ggplot(primary, aes(x = observed_day, y = predicted_day, colour = group_simple)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#888888", linewidth = 0.75) +
  geom_point(size = 3.0, alpha = 0.90,
             position = position_jitter(width = 0.08, height = 0.08, seed = 20260824)) +
  scale_colour_manual(values = palette, drop = FALSE) +
  scale_x_continuous(breaks = c(3, 5, 7, 9, 11), limits = c(2.8, 11.2)) +
  scale_y_continuous(breaks = c(3, 5, 7, 9, 11), limits = c(2.8, 11.2)) +
  annotate("text", x = 10.95, y = 4.13, hjust = 1, vjust = 0,
           label = paste0(
             "rho = ", sprintf("%.2f", rho), "\n",
             "MAE = ", sprintf("%.2f", mae), " days\n",
             "Exact stage = ", sprintf("%.0f%%", 100 * exact), "\n",
             "Pairwise order = ", sprintf("%.1f%%", 100 * pairwise)
           ), size = 4.0, lineheight = 1.08, colour = "#202020") +
  labs(x = "Observed luteal day", y = "Predicted luteal day", colour = "group_simple") +
  theme_figure1c

save_plot <- function(plot, base, width, height) {
  ggsave(paste0(base, ".png"), plot = plot, width = width, height = height,
         units = "in", dpi = 600, bg = "white", limitsize = FALSE)
  ggsave(paste0(base, ".tiff"), plot = plot, width = width, height = height,
         units = "in", dpi = 600, bg = "white", limitsize = FALSE,
         device = "tiff", compression = "lzw")
  ggsave(paste0(base, ".pdf"), plot = plot, width = width, height = height,
         units = "in", bg = "white", limitsize = FALSE, device = cairo_pdf)
}

# The current locked standalone Figure 1C is a single fertile validation plot.
# No panel letter is added here; the output is intended for author-controlled replacement.
save_plot(figure1c, file.path(output_root, "01_Figure1C", "Figure1_panelC_corrected"),
          width = 3750 / 600, height = 3298 / 600)

# Supplementary S1: preserve the current clean standalone geometry and labels.
s1 <- primary
s1$group_simple <- factor(as.character(s1$group_simple),
                           levels = c("Fertile_LH7", "LH11", "LH3", "LH5", "LH9"))
s1_plot <- ggplot(s1, aes(x = observed_day, y = predicted_day, colour = group_simple)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#888888", linewidth = 0.75) +
  geom_point(size = 1.85, alpha = 0.90,
             position = position_jitter(width = 0.08, height = 0.08, seed = 20260824)) +
  scale_colour_manual(values = palette, drop = FALSE) +
  scale_x_continuous(breaks = c(3, 5, 7, 9, 11), limits = c(2.8, 11.2)) +
  scale_y_continuous(breaks = c(3, 5, 7, 9, 11), limits = c(2.8, 11.2)) +
  annotate("text", x = 10.82, y = 4.13, hjust = 1, vjust = 0,
           label = paste0(
             "LOOCV rho = ", sprintf("%.2f", rho), "\n",
             "MAE = ", sprintf("%.2f", mae), " days\n",
             "Exact stage = ", sprintf("%.0f%%", 100 * exact)
           ), size = 2.35, lineheight = 1.05, colour = "#202020") +
  labs(x = "Observed luteal day", y = "Predicted luteal day", colour = "group_simple") +
  theme_bw(base_size = 7.8, base_family = "Helvetica") +
  theme(
    panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.30),
    panel.grid.minor = element_line(colour = "#F0F0F0", linewidth = 0.20),
    panel.border = element_rect(colour = "#555555", linewidth = 0.45, fill = NA),
    axis.line = element_line(colour = "#303030", linewidth = 0.35),
    axis.ticks = element_line(colour = "#303030", linewidth = 0.35),
    axis.text = element_text(size = 6.4),
    axis.title = element_text(size = 7.8),
    legend.title = element_text(size = 7.2),
    legend.text = element_text(size = 6.8),
    legend.key = element_blank(),
    plot.margin = margin(5, 10, 5, 7, unit = "pt")
  )
save_plot(s1_plot, file.path(output_root, "02_Supplementary_S1", "Supplementary_Figure_S1_corrected"),
          width = 1924 / 600, height = 1508 / 600)

# Supplementary S3: preserve historical violin/boxplot structure and dimensions,
# replacing the historical first label with the corrected data-derived control.
low_assoc_null$control <- "Data-derived low-association gene-set control"
permuted_null$control <- "Permuted LH labels"
random_null$control <- "Random 80-gene sets"
null_runs <- rbind(
  low_assoc_null[, c("control", "rho")],
  permuted_null[, c("control", "rho")],
  random_null[, c("control", "rho")]
)
null_runs$control <- factor(null_runs$control,
                            levels = c("Data-derived low-association gene-set control",
                                       "Permuted LH labels", "Random 80-gene sets"))
s3_palette <- c(
  "Data-derived low-association gene-set control" = "#E0E0E0",
  "Permuted LH labels" = "#F4A582",
  "Random 80-gene sets" = "#92C5DE"
)
s3_labels <- c(
  "Data-derived low-association gene-set control" = "Data-derived low-association\ngene-set control",
  "Permuted LH labels" = "Permuted\nLH labels",
  "Random 80-gene sets" = "Random 80-gene\nsets"
)
s3_plot <- ggplot(null_runs, aes(x = control, y = rho, fill = control)) +
  geom_violin(alpha = 0.30, colour = NA, width = 0.92, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, colour = "grey25",
               fill = "white", linewidth = 0.55) +
  geom_hline(yintercept = rho, linetype = 2, colour = "#D73027", linewidth = 0.8) +
  annotate("text", x = 1.5, y = 0.82, hjust = 0.5, vjust = 0,
           label = paste0("Observed 80-gene\nLOOCV rho = ", sprintf("%.2f", rho)),
           colour = "#BD1F36", size = 2.5, lineheight = 0.95) +
  scale_fill_manual(values = s3_palette, drop = FALSE) +
  scale_x_discrete(labels = s3_labels) +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
  labs(x = NULL, y = "LOOCV Spearman rho") +
  theme_bw(base_size = 7.8, base_family = "Helvetica") +
  theme(
    panel.grid.major = element_line(colour = "#E6E6E6", linewidth = 0.30),
    panel.grid.minor = element_line(colour = "#F0F0F0", linewidth = 0.20),
    panel.border = element_rect(colour = "#555555", linewidth = 0.45, fill = NA),
    axis.line = element_line(colour = "#303030", linewidth = 0.35),
    axis.ticks = element_line(colour = "#303030", linewidth = 0.35),
    axis.text.x = element_text(size = 5.8, colour = "#303030", lineheight = 0.9),
    axis.text.y = element_text(size = 6.5, colour = "#303030"),
    axis.title.y = element_text(size = 7.5, colour = "#202020"),
    legend.position = "none",
    plot.margin = margin(5, 7, 6, 7, unit = "pt")
  )
save_plot(s3_plot, file.path(output_root, "03_Supplementary_S3", "Supplementary_Figure_S3_corrected"),
          width = 2288 / 600, height = 1466 / 600)

writeLines(c(
  "CORRECTED FIGURE REBUILD REPORT",
  "",
  "Outputs:",
  paste0("Figure 1C corrected panel: ", file.path(output_root, "01_Figure1C", "Figure1_panelC_corrected.{png,tiff,pdf}")),
  paste0("Supplementary S1 corrected: ", file.path(output_root, "02_Supplementary_S1", "Supplementary_Figure_S1_corrected.{png,tiff,pdf}")),
  paste0("Supplementary S3 corrected: ", file.path(output_root, "03_Supplementary_S3", "Supplementary_Figure_S3_corrected.{png,tiff,pdf}")),
  "",
  "Corrected input files:",
  pred_path,
  summary_path,
  fold_audit_path,
  random_path,
  permuted_path,
  low_assoc_path,
  null_summary_path,
  "",
  "Historical style baselines inspected:",
  "outputs/Figure_1C_corrected_nested_LOOCV.png",
  "outputs/Figure_S1_corrected_nested_LOOCV.png",
  "outputs/Figure_S3_corrected_negative_controls.png",
  "scripts/internal_validation/rebuild_corrected_figure1c_s1_s3_public.R",
  "scripts/internal_validation/corrected_nested_validation_and_nulls_public.R",
  "The corrected panel style was reproduced from the locked standalone assets; no historical figure was overwritten.",
  "",
  "Figure 1C panel letter: none added (standalone output; no C label).",
  "Figure 1 composite preview: NOT GENERATED per task restriction.",
  "",
  "Corrected primary metrics used:",
  paste0("rho = ", format(rho, digits = 16)),
  paste0("MAE = ", format(mae, digits = 16)),
  paste0("pairwise concordance = ", format(pairwise, digits = 16)),
  paste0("exact nearest-stage agreement = ", format(exact, digits = 16)),
  "",
  "Corrected null-control inputs:",
  "Random 80-gene sets: corrected_random80_null_metrics.csv",
  "Permuted LH labels: corrected_permutedLH_null_metrics.csv",
  "Data-derived low-association gene-set control: corrected_low_association_null_metrics.csv",
  "Random median rho = 0.753852; 95% null interval = 0.614612 to 0.864465; max rho = 0.884266; b = 0/250; plus-one P = 0.003984.",
  "Permuted LH median rho = -0.122461; 95% null interval = -0.623042 to 0.376688; max rho = 0.537557; b = 0/250; plus-one P = 0.003984.",
  "Data-derived low-association median rho = 0.432590; 95% null interval = 0.034989 to 0.689521; max rho = 0.769757; b = 0/250; plus-one P = 0.003984.",
  "",
  "Formal files modified: NO.",
  "Historical figures overwritten: NO.",
  "No manuscript, supplementary appendix, tables, legends, formal submission files, GitHub, or Zenodo files modified.",
  "No statistical analysis rerun; corrected result files were read only."
), file.path(output_root, "CORRECTED_FIGURE_REBUILD_REPORT.md"))

message("CORRECTED_FIGURE_ASSETS_COMPLETE")
