#!/usr/bin/env python3
"""Index verified disease-level meta-analysis exports; no recalculation."""
from pathlib import Path
import csv

BASE = Path(__file__).resolve().parents[1]
for name in ["meta_rif_verified.tsv", "meta_endometriosis_verified.tsv", "meta_adenomyosis_verified.tsv"]:
    with (BASE / "source_data" / "verified_sensitivity" / name).open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert rows
    print(name, len(rows))
