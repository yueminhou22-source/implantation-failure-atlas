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

rows = read_csv("source_data/external_validation/GSE144895_timing_scores.csv")
preg = [float(r["timing_score"]) for r in rows if r["pregnancy_group"] == "Pregnant"]
non = [float(r["timing_score"]) for r in rows if r["pregnancy_group"] == "Non-pregnant"]
delta = statistics.mean(preg) - statistics.mean(non)

# Exact two-sided Mann-Whitney U P value by enumeration.
# There are C(20,9)=167,960 assignments, which is small enough for direct enumeration.
all_scores = preg + non
ranks = rankdata(all_scores)
m, n = len(preg), len(non)
obs_rank_sum = sum(ranks[:m])
u_obs = obs_rank_sum - m*(m+1)/2

# The released scores are distinct, so the null U distribution depends only on choosing m ranks from 1..m+n.
N = m + n
counts = {}
total = 0
for comb in itertools.combinations(range(1, N+1), m):
    u = sum(comb) - m*(m+1)/2
    counts[u] = counts.get(u, 0) + 1
    total += 1

lower = sum(c for u,c in counts.items() if u <= u_obs) / total
upper = sum(c for u,c in counts.items() if u >= u_obs) / total
p_exact = min(1.0, 2 * min(lower, upper))

print(f"patients={len(rows)}")
print(f"pregnant_n={m}")
print(f"non_pregnant_n={n}")
print(f"pregnant_minus_nonpregnant_mean_delta={delta:.12f}")
print(f"mann_whitney_u={u_obs:.6f}")
print(f"exact_two_sided_mann_whitney_p={p_exact:.15g}")
