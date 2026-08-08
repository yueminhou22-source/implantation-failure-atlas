#!/usr/bin/env python3
"""Check the secondary rank-based timing validation export; no recalculation."""
from pathlib import Path
import csv

BASE = Path(__file__).resolve().parents[1]
with (BASE / "source_data" / "verified_sensitivity" / "rank_based_evidence_summary.tsv").open(encoding="utf-8", newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))
assert rows
print(f"rank-based evidence rows: {len(rows)}")
