# External-validation reconstruction scripts

These four standalone Python scripts reconstruct the primary point/rank statistics from the released timing-score CSV files:

- `reproduce_GSE234368_summary.py`
- `reproduce_GSE98386_summary.py`
- `reproduce_GSE180485_summary.py`
- `reproduce_GSE144895_summary.py`

They use only the Python standard library.

Run from any working directory, for example:

```bash
python3 scripts/external_validation/reproduce_GSE234368_summary.py
python3 scripts/external_validation/reproduce_GSE98386_summary.py
python3 scripts/external_validation/reproduce_GSE180485_summary.py
python3 scripts/external_validation/reproduce_GSE144895_summary.py
```

These scripts operate on the released row-level timing-score tables. They do not reconstruct the upstream raw-expression preprocessing or the frozen-score derivation from raw GEO files, and they do not re-bootstrap the confidence intervals preserved in Supplementary Table S15.
