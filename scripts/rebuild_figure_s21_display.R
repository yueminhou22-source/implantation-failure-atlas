suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

root <- normalizePath(".", mustWork = TRUE)
source_path <- file.path(root, "source_data", "heca_non_endo_key_shift_comparison.csv")
out_dir <- file.path(root, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

source_table <- read_csv(source_path, show_col_types = FALSE)

# Replace internal state labels only for display; estimates and ordering remain unchanged.
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

source_table <- source_table %>% mutate(display_celltype = display_celltype(celltype))

plot_data <- bind_rows(
  source_table %>% transmute(
    label = paste(dataset, display_celltype, sep = " | "),
    reference = "Full HECA",
    estimate = estimate_full,
    lower = lower_full,
    upper = upper_full
  ),
  source_table %>% transmute(
    label = paste(dataset, display_celltype, sep = " | "),
    reference = "Non-endometriosis HECA",
    estimate = estimate_non_endo,
    lower = lower_non_endo,
    upper = upper_non_endo
  )
) %>%
  mutate(label = factor(label, levels = rev(unique(paste(source_table$dataset, source_table$display_celltype, sep = " | ")))))

plot <- ggplot(plot_data, aes(estimate, label, color = reference)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey70") +
  geom_errorbar(aes(xmin = lower, xmax = upper), orientation = "y", position = position_dodge(width = 0.55), width = 0.18, linewidth = 0.75) +
  geom_point(position = position_dodge(width = 0.55), size = 3) +
  scale_color_manual(values = c("Full HECA" = "#F8766D", "Non-endometriosis HECA" = "#00BFC4")) +
  labs(
    title = "Key cell-state shifts: full versus non-endometriosis HECA",
    x = "Predicted-day shift",
    y = NULL,
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", margin = margin(b = 8)),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(12, 28, 12, 12)
  )

ggsave(file.path(out_dir, "Figure_S21_rebuilt_display_labels.png"), plot, width = 8.6, height = 5.1, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Figure_S21_rebuilt_display_labels.tiff"), plot, width = 8.6, height = 5.1, dpi = 600, compression = "lzw", limitsize = FALSE)
ggsave(file.path(out_dir, "Figure_S21_rebuilt_display_labels.pdf"), plot, width = 8.6, height = 5.1, limitsize = FALSE)