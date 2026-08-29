# Release notes — v1.2.0

## Corrected nested internal validation

This release adds the corrected, isolated nested fertile-sample leave-one-out validation and its three prespecified null controls. The corrected implementation uses training-only filtering, association testing, BH adjustment, feature selection, orientation, centering/scaling and calibration for each held-out fold.

The corrected primary summary is Spearman rho = 0.925616, MAE = 0.956520 days, pairwise concordance = 0.960317 and exact nearest-stage agreement = 0.500000 across 18/18 estimable held-out samples.

Random 80-gene sets, permuted LH labels and the data-derived low-association gene-set control each contain B = 250 estimable iterations. The plus-one empirical P value is 0.003984 for each control in the corrected summary. The data-derived low-association control is the current label for the historical housekeeping-like control; no canonical housekeeping list is asserted by this release.

## Scope

The frozen 80-gene analytical timing coordinate and the independent external evaluation materials added in v1.1.0 are not redefined by this validation update. Corrected Figure 1C, S1 and S3 assets are included as source-table-derived outputs. No raw expression matrices are mirrored.
