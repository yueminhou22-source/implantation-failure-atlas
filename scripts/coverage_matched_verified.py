#!/usr/bin/env python3
"""Check the verified original-versus-coverage-matched export; no recalculation."""
from pathlib import Path
import csv

BASE = Path(__file__).resolve().parents[1]
path = BASE / "source_data" / "verified_sensitivity" / "verified_original_vs_coverage_matched.tsv"
with path.open(encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
assert rows
assert all(abs(float(r["score_difference"])) < 1e-10 for r in rows if r.get("score_difference") not in (None, ""))
print(f"verified rows: {len(rows)}")
