#!/usr/bin/env python3
"""Lightweight reconstruction from released row-level timing-score source data."""
from pathlib import Path
import csv, math, statistics, itertools

ROOT = Path(__file__).resolve().parents[2]

def read_csv(relpath):
    with (ROOT / relpath).open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))

def rankdata(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    k = 0
    while k < len(order):
        j = k + 1
        while j < len(order) and values[order[j]] == values[order[k]]:
            j += 1
        avg = (k + 1 + j) / 2.0
        for p in range(k, j):
            ranks[order[p]] = avg
        k = j
    return ranks

def pearson(x, y):
    mx, my = statistics.mean(x), statistics.mean(y)
    num = sum((a-mx)*(b-my) for a,b in zip(x,y))
    den = math.sqrt(sum((a-mx)**2 for a in x) * sum((b-my)**2 for b in y))
    return num / den

def spearman(x, y):
    return pearson(rankdata(x), rankdata(y))

rows = read_csv("source_data/external_validation/GSE180485_timing_scores.csv")
included = [r for r in rows if r["included_in_correlation"].strip().upper() == "TRUE" and r["reported_LH_day"].strip() != ""]
lh = [float(r["reported_LH_day"]) for r in included]
scores = [float(r["timing_score"]) for r in included]
rho = spearman(lh, scores)

print(f"independent_biopsies={len(rows)}")
print(f"correlation_n={len(included)}")
print(f"spearman_rho={rho:.12f}")
