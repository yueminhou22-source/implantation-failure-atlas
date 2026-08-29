# Figure S3 v1.2.0 corrected asset note

The current Figure S3 source mapping uses the corrected nested-validation null-control tables under `source_data/internal_validation/` and the public rebuild script under `scripts/internal_validation/rebuild_corrected_figure1c_s1_s3_public.R`.

The three current control labels are:

1. Data-derived low-association gene-set control;
2. Permuted LH labels;
3. Random 80-gene sets.

The historical housekeeping-like control and the historical LOOCV metrics are not current v1.2.0 validation outputs. The corrected source tables preserve all 250 iteration rows for each control and the corrected observed metrics. The corrected figure assets are source-table-derived outputs and do not modify the formal manuscript or supplementary appendix.
