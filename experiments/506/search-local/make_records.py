"""Collect the best configuration per n (n = 9..16) from runs/sa_results.jsonl, certify it exactly
with verify_exact.certify, and write runs/records.json.  Rational (integer-scalable) records are
preferred when several configurations tie."""
import json, sys, traceback
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/search-local')
from verify_exact import certify
HERE = '/home/user/Solving-hard-math-problems-/experiments/506/search-local'


def formula(n): return (n - 1) * (n - 2) // 2 + 1 - (n - 1) // 2


res = []
for line in open(f"{HERE}/runs/sa_results_checked.jsonl") if __import__('os').path.exists(f"{HERE}/runs/sa_results_checked.jsonl") else open(f"{HERE}/runs/sa_results.jsonl"):
    try:
        res.append(json.loads(line))
    except Exception:
        pass
by_n = {}
for r in res:
    by_n.setdefault(r['n'], []).append(r)
records = {}
for n in sorted(by_n):
    cands = sorted(by_n[n], key=lambda r: (r['true'], r['exact'] is None, not all(isinstance(c, list) and all('/' in x or x.lstrip('-').isdigit() for x in c) for c in (r['exact'] or []) if c != 'inf')))
    best_val = cands[0]['true']
    out = None
    for r in cands:
        if r['true'] != best_val or r['exact'] is None:
            continue
        try:
            cert = certify(r['exact'], r['centre'], want=r['true'], verbose=False)
        except Exception as e:
            print(f"n={n} {r['universe']}: certification failed: {e!r}")
            continue
        if cert['count'] != r['true']:
            print(f"n={n} {r['universe']}: certified count {cert['count']} != claimed {r['true']}")
            continue
        out = dict(n=n, circles=cert['count'], formula=formula(n), universe=r['universe'], structure=r['structure'],
                   certified=cert, moebius_S=r['exact'], centre=r['centre'])
        print(f"n={n:2d} f={formula(n):3d} best={cert['count']:3d} universe={r['universe']} exact={cert['exact']} coords={cert['coords']}")
        break
    if out is None:
        print(f"n={n}: best {best_val} could not be certified from stored exact data")
        out = dict(n=n, circles=best_val, formula=formula(n), universe=cands[0]['universe'], structure=cands[0]['structure'], certified=None)
    records[n] = out
json.dump(records, open(f"{HERE}/runs/records.json", "w"), indent=1, default=str)
