#!/usr/bin/env python3
"""Validate the source-table reproducibility package without external data."""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
NEEDED = [
    "source_data/fertile_day_spearman_round5.csv",
    "source_data/receptivity_signature_genes_round5.csv",
    "source_data/internal_validation/corrected_primary_loocv_predictions.csv",
    "source_data/internal_validation/corrected_primary_loocv_summary.csv",
    "source_data/internal_validation/corrected_random80_null_metrics.csv",
    "source_data/internal_validation/corrected_permutedLH_null_metrics.csv",
    "source_data/internal_validation/corrected_low_association_null_metrics.csv",
    "source_data/internal_validation/corrected_negative_control_summary.csv",
    "source_data/heca_mapping_quality_heatmap_values.csv",
]
missing = [item for item in NEEDED if not (ROOT / item).exists()]
if missing:
    print("Missing required source-table inputs:", *missing, sep="\n- ")
    sys.exit(1)
print("Input check passed. Source-table reproduction scripts can be run from:", ROOT)
