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

rows = read_csv("source_data/external_validation/GSE98386_paired_timing_scores.csv")
d = [float(r["paired_delta"]) for r in rows]
nonzero = [x for x in d if x != 0]
positive = sum(x > 0 for x in nonzero)

# Exact two-sided Wilcoxon signed-rank P for the released data.
# The 20 non-zero paired differences are all positive, so the observed
# signed-rank sum is at the extreme tail: p = 2 / 2^n.
if positive == len(nonzero):
    p_exact = 2.0 / (2 ** len(nonzero))
else:
    raise RuntimeError("Released paired deltas are no longer all positive; use a full exact signed-rank implementation.")

print(f"women={len(rows)}")
print(f"mean_paired_delta={statistics.mean(d):.12f}")
print(f"positive_pairs={positive}/{len(nonzero)}")
print(f"exact_two_sided_wilcoxon_p={p_exact:.15g}")
