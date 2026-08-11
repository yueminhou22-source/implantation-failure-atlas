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

rows = read_csv("source_data/external_validation/GSE234368_timing_scores.csv")
stages = [int(r["reported_stage"]) for r in rows]
scores = [float(r["timing_score"]) for r in rows]
rho = spearman(stages, scores)

ordered = 0.0
total = 0
for i in range(len(rows)):
    for j in range(i+1, len(rows)):
        si, sj = stages[i], stages[j]
        if si == sj:
            continue
        total += 1
        low_score, high_score = (scores[i], scores[j]) if si < sj else (scores[j], scores[i])
        if low_score < high_score:
            ordered += 1.0
        elif low_score == high_score:
            ordered += 0.5

counts = {s: stages.count(s) for s in sorted(set(stages))}
medians = {s: statistics.median([scores[i] for i,v in enumerate(stages) if v == s]) for s in counts}

print(f"n={len(rows)}")
print("stage_counts=" + ", ".join(f"{k}:{v}" for k,v in counts.items()))
print(f"spearman_rho={rho:.12f}")
print(f"ordered_pair_concordance_percent={100*ordered/total:.10f}")
print("stage_medians=" + ", ".join(f"{k}:{v:.12f}" for k,v in medians.items()))
