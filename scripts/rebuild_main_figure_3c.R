#!/usr/bin/env Rscript

# Rebuilds the Figure 3C correlation heatmap from the included source table.
suppressPackageStartupMessages({ library(ggplot2); library(readr) })
all_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", all_args[grep("^--file=", all_args)][1])
script_dir <- dirname(normalizePath(script_path))
root <- normalizePath(file.path(script_dir, ".."))
src <- file.path(root, "source_data", "heca_mapping_quality_heatmap_values.csv")
out <- file.path(root, "outputs")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(src)) stop("Missing source table: ", src)
df <- read_csv(src, show_col_types = FALSE)
df$dataset <- factor(df$dataset, levels = c("GSE250130", "GSE179640", "GSE213216", "GSE214411"))
df$celltype <- factor(df$celltype, levels = c("Luminal_Epi", "Glandular_Epi", "Stroma", "Decidual_Stroma", "Endothelial", "Lymphoid", "Myeloid"))
p <- ggplot(df, aes(celltype, dataset, fill = mean_cor)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", mean_cor)), size = 3.0) +
  scale_fill_gradient(low = "#E6F0FA", high = "#144B7D", name = "Mean r") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 10, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 26, hjust = 1), plot.margin = margin(7, 12, 7, 7, "pt"))
ggsave(file.path(out, "Figure_3C_rebuilt_from_source_table.png"), p, width = 5.4, height = 4.2, dpi = 600, limitsize = FALSE)
ggsave(file.path(out, "Figure_3C_rebuilt_from_source_table.pdf"), p, width = 5.4, height = 4.2, limitsize = FALSE)