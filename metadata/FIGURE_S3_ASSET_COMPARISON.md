# Supplementary Figure S3 Asset Comparison

## Scope

This audit compared the legacy FINAL_DATA_AUDITED assets with the authoritative S3 rebuild generated from the verified negative-control source table and the current source-table rebuild script.

## Legacy discrepancy

| Asset | SHA-256 | Dimensions | DPI | Status |
|---|---|---:|---:|---|
| Legacy DOCX embedded Figure S3 | `d6ed7a34...` | 5640 x 3360 px | 600 | Not synchronized with the legacy figure-package PNG/TIFF. |
| Legacy Additional file 2 PNG | `5b82903e...` | 2288 x 1390 px | not embedded | Different raster asset. |
| Legacy Additional file 2 TIFF | `780ebf9d...` | 2288 x 1390 px | 300 | Different raster asset. |

## Authoritative rebuilt asset

The authoritative Figure S3 was rebuilt from `source_data/fertile_timeline_negative_controls.csv` and `source_data/fertile_timeline_negative_control_summary_verified.csv` by `scripts/rebuild_audit_supplementary_figures.R`.

| Final asset | SHA-256 | Dimensions | DPI | Status |
|---|---|---:|---:|---|
| Figure_S3.png | `c6581be73e76e5c37b4293820522a7753d0d0e6ad3898a828830c115175a9b9d` | 5640 x 3360 px | 600 | Used in Additional file 2 and embedded in the final supplementary appendix. |
| Figure_S3.tiff | `ef69bae65249787d038da3210264b4a71bcaf96f51ecc7e1693895d758577d31` | 5640 x 3360 px | 600 | Used in Additional file 2. |

## Visual and content checks

- The rebuilt figure retains the three required controls: housekeeping-like gene sets, permuted LH labels, and random 80-gene sets.
- The observed LOOCV-rho reference line is retained and fully visible.
- The figure contains no long embedded caption; explanatory text remains in the Word caption.
- No top, right, or bottom clipping was observed in the rendered supplementary appendix.
- The figure legend continues to report the verified null result: 0/250 null rho values greater than or equal to the observed rho and plus-one empirical P = 0.004.

## Conclusion

The final DOCX embedded PNG and the Additional file 2 PNG are the same authoritative rebuild. The final TIFF was generated in the same rebuild run from the same plotting object and source tables.