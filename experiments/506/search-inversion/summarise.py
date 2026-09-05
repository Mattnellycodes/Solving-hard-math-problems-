"""Summary of runs/results.jsonl and runs/records.jsonl: best count per n, and the Moebius block
structure (block sizes, largest block, whether an (n-1)-block is present) of every record with
count <= f(n)."""
import json, os, sys
from collections import defaultdict
import numpy as np
import sphere as SP

HERE = os.path.dirname(os.path.abspath(__file__))


def load(path):
    out = []
    if os.path.exists(path):
        for line in open(path):
            try:
                out.append(json.loads(line))
            except Exception:
                pass
    return out


def structure_of_record(rec):
    V = np.array(rec['S_xyz'], float)
    ev = SP.evaluate(V)
    if ev is None:
        return None
    o = np.array(rec['best_centre_xyz'], float)
    thr = SP.blocks_through(o, ev.Nn, ev.D)
    sizes = defaultdict(int)
    for b in ev.blocks:
        sizes[len(b)] += 1
    lines = sorted(len(ev.blocks[t]) for t in thr)
    return dict(sizes=dict(sorted(sizes.items())), nblocks=ev.nb, deg=len(thr), line_sizes=lines,
                largest=max(sizes), has_big=(max(sizes) >= len(V) - 1), count=ev.nb - len(thr))


if __name__ == "__main__":
    res = load(os.path.join(HERE, 'runs', 'results.jsonl'))
    recs = load(os.path.join(HERE, 'runs', 'records.jsonl'))
    best = {}
    per_univ = defaultdict(dict)
    for r in res:
        n = r['n']; per_univ[r['universe']][n] = r['count']
        if n not in best or r['count'] < best[n]['count']:
            best[n] = r
    print("=== best per n (all universes so far) ===")
    for n in sorted(best):
        r = best[n]
        flag = "BELOW" if r['count'] < r['formula'] else ("=f" if r['count'] == r['formula'] else f"+{r['count']-r['formula']}")
        print(f"n={n:2d} f={r['formula']:3d} best={r['count']:3d} [{flag}] {r['universe']} centre={r['centre_label']}")
    print("\n=== per universe (count - f(n)) ===")
    for u in sorted(per_univ):
        d = per_univ[u]
        print(f"{u:28s} " + " ".join(f"n{n}:{d[n]-SP.formula(n):+d}" for n in sorted(d)))
    if '--records' in sys.argv:
        print("\n=== records (count <= f(n)) with structure ===")
        seen = set()
        for rec in recs:
            key = (rec['universe'], rec['n'])
            if key in seen:
                continue
            seen.add(key)
            s = structure_of_record(rec)
            print(f"{rec['universe']:22s} n={rec['n']:2d} count={rec['count']} f={rec['formula']} structure={s}")
