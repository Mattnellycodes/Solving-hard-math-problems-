"""Best count per n over all SA runs (runs/sa_results.jsonl) and per-universe table."""
import json, sys
from collections import defaultdict
HERE = '/home/user/Solving-hard-math-problems-/experiments/506/search-local'
def formula(n): return (n - 1) * (n - 2) // 2 + 1 - (n - 1) // 2
res = []
for line in open(f"{HERE}/runs/sa_results_checked.jsonl") if __import__('os').path.exists(f"{HERE}/runs/sa_results_checked.jsonl") else open(f"{HERE}/runs/sa_results.jsonl"):
    try: res.append(json.loads(line))
    except Exception: pass
best = {}
per = defaultdict(dict)
for r in res:
    n = r['n']
    if n not in best or r['true'] < best[n]['true']:
        best[n] = r
    u = r['universe']
    if n not in per[u] or r['true'] < per[u][n]:
        per[u][n] = r['true']
print("=== best per n ===")
for n in sorted(best):
    r = best[n]
    flag = "BELOW" if r['true'] < formula(n) else ("=f" if r['true'] == formula(n) else f"+{r['true']-formula(n)}")
    print(f"n={n:2d} f={formula(n):3d} best={r['true']:3d} [{flag}] {r['universe']} structure={r['structure']}")
print("\n=== per universe: count - f(n) ===")
for u in sorted(per):
    d = per[u]
    print(f"{u:18s} " + " ".join(f"n{n}:{d[n]-formula(n):+d}" for n in sorted(d)))
