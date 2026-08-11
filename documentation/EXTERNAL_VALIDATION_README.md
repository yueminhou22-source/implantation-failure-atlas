# Independent external validation of the frozen analytical timing coordinate

This directory documents the independent external-validation materials added for release v1.1.0.

## Frozen timing coordinate

The analytical timing coordinate was defined before evaluation of these external cohorts as:

`timing score = mean(late genes) - mean(early genes)`

using the prespecified 40 late and 40 early genes in `source_data/receptivity_signature_genes_round5.csv`.

For each external cohort, unavailable signature genes were omitted without replacement or imputation. No gene reselection, model refitting or target-cohort recalibration was performed.

The row-level files in `source_data/external_validation/` contain the final cohort-level timing scores used for the reported external-validation analyses. Raw GEO expression matrices are not mirrored in this repository; the original public repositories remain the authoritative raw-data source.

## Cohorts

- **GSE234368** — 150 samples; Stage 5 n=51, Stage 6 n=66, Stage 7 n=33; 54/80 frozen genes available (24 early, 30 late). Primary support: Spearman ordering across Stage 5–7 and ordered-pair concordance.
- **GSE98386** — 20 women with paired LH+2 and LH+8 biopsies; 68/80 genes available (30 early, 38 late). Primary support: within-person LH+8-minus-LH+2 timing-score change.
- **GSE180485** — 36 independent biopsies; 63/80 genes available (29 early, 34 late). Thirty-five biopsies had source-reported LH+ timing and were included in the correlation analysis.
- **GSE144895** — 20 patients (9 pregnant, 11 non-pregnant); 68/80 genes available (30 early, 38 late). This cohort was retained as an exploratory outcome-linked assessment and was not used to establish clinical prediction or utility.

## Released files

- `GSE234368_timing_scores.csv`
- `GSE98386_paired_timing_scores.csv`
- `GSE180485_timing_scores.csv`
- `GSE144895_timing_scores.csv`
- `Supplementary_Table_S15.csv`
- `frozen_gene_coverage.csv`
- `external_validation_manifest.csv`

## Lightweight reconstruction scripts

The scripts in `scripts/external_validation/` reconstruct the primary point/rank statistics directly from the released row-level timing-score tables. They are intentionally scoped to the released source tables and do not constitute a raw-GEO-to-score processing pipeline.

The confidence intervals reported in Supplementary Table S15 are preserved as frozen analysis outputs in `Supplementary_Table_S15.csv`; the lightweight scripts do not re-bootstrap those intervals.

## Interpretation boundary

These materials support transferability of the analytical timing coordinate across independent datasets. They do not establish absolute physiological dating, a clinically deployable window-of-implantation test, clinical prediction, or clinical utility.
