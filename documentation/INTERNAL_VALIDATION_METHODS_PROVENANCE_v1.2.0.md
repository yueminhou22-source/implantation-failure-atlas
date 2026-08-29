# Internal validation methods and provenance — v1.2.0

This note documents the corrected nested validation update included in v1.2.0. It is written with repository-relative paths; raw expression matrices and local execution paths are not mirrored in this repository.

## Scope and frozen-coordinate boundary

The full-cohort frozen 80-gene analytical timing coordinate is unchanged. This update evaluates its fertile-sample predictive ordering under a nested leave-one-sample-out procedure and does not redefine the frozen gene membership, the historical full-cohort coordinate or the independent external evaluation materials added in v1.1.0.

## Corrected primary procedure

For each of the 18 fertile samples, all filtering, gene–day association, Benjamini–Hochberg adjustment, feature selection, early/late orientation, centering/scaling and linear calibration were performed using the 17 training samples only. The score was defined as `mean(late) − mean(early)`, and the training-only calibration was `day ~ score`. Training-derived filtering and scaling parameters were then applied to the held-out sample. The locked fold audit is provided in `source_data/internal_validation/corrected_primary_loocv_fold_audit.csv` and the sample-level predictions are provided in `source_data/internal_validation/corrected_primary_loocv_predictions.csv`.

Corrected summary: Spearman rho = 0.925616344099057; MAE = 0.956520006069255 days; pairwise concordance = 0.960317460317460; exact nearest-stage agreement = 0.500000; 18/18 estimable folds.

## Negative controls

The release provides 250 iterations for each of:

- `random_80_gene`;
- `permuted_LH`;
- `data_derived_low_association`.

The current public name for the third control is “data-derived low-association gene-set control”; it is not a canonical housekeeping list. Each control has line-level iteration metrics in `source_data/internal_validation/` and the sampled-gene/fold provenance is in `source_data/internal_validation/corrected_negative_control_fold_provenance.csv`. The corrected summary reports the median rho, 2.5th–97.5th percentile interval, B, exceedance count and plus-one empirical P.

The corrected summaries are:

| Control | Median rho | 2.5th–97.5th interval | B | Exceedance count | Plus-one P |
|---|---:|---:|---:|---:|---:|
| Random 80-gene sets | 0.753852 | 0.614612 to 0.864465 | 250 | 0 | 0.003984 |
| Permuted LH labels | -0.122461 | -0.623042 to 0.376688 | 250 | 0 | 0.003984 |
| Data-derived low-association gene-set control | 0.432590 | 0.034989 to 0.689521 | 250 | 0 | 0.003984 |

The plus-one empirical P is `(b + 1) / (B + 1)`, implemented directly in the corrected summary path.

## Public files and scripts

- Primary prediction table: `source_data/internal_validation/corrected_primary_loocv_predictions.csv`.
- Primary summary: `source_data/internal_validation/corrected_primary_loocv_summary.csv`.
- Fold audit: `source_data/internal_validation/corrected_primary_loocv_fold_audit.csv`.
- Null metrics: `source_data/internal_validation/corrected_random80_null_metrics.csv`, `corrected_permutedLH_null_metrics.csv`, and `corrected_low_association_null_metrics.csv`.
- Null summary: `source_data/internal_validation/corrected_negative_control_summary.csv`.
- Null fold provenance: `source_data/internal_validation/corrected_negative_control_fold_provenance.csv`.
- Validation execution source: `scripts/internal_validation/corrected_nested_validation_and_nulls_public.R`.
- Corrected figure rebuild source: `scripts/internal_validation/rebuild_corrected_figure1c_s1_s3_public.R`.

The public execution source expects raw/input paths to be supplied by the user at run time and does not embed local author paths. Runtime/package details are summarized in the public environment files.
