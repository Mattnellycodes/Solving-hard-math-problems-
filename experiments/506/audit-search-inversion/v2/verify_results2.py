"""Exact re-verification of every S2lat row of search-inversion/runs/results2.jsonl whose reported count
equals f(n) (n >= 10).  For each row the point set S is rebuilt from S_labels; the best projection centre is
found exactly among the sphere points on >= 2 block planes of S (rational ones; irrational conjugate
pairs handled as lines); the circle count R - L is cross-checked by exact inversion + coplanar count."""
import json, re, sys
from fractions import Fraction as Fr
from itertools import combinations
from math import comb
from collections import defaultdict, Counter
from s2lat_check import invert, count_coplanar, cross, sub, dot, prim
from sphsearch import plane_key, line_of, sqrt_fr, contains, f

def best_centre(S, N):
    planes = defaultdict(set)
    for i, j, k in combinations(range(len(S)), 3):
        key = plane_key(S[i], S[j], S[k])
        if key: planes[key].update((i, j, k))
    plist = list(planes); R = len(plist)
    best = (1, None)
    for k1, k2 in combinations(plist, 2):
        ln = line_of(k1, k2)
        if ln is None: continue
        p0, u = ln
        A = Fr(dot(u, u)); B = 2 * dot(p0, u); C = dot(p0, p0) - N
        disc = B * B - 4 * A * C
        if disc < 0: continue
        s = sqrt_fr(disc)
        cs = [('irr', p0, u)] if s is None else [('rat', tuple(p + ((-B + sg * s) / (2 * A)) * x for p, x in zip(p0, u))) for sg in (1, -1)]
        for c in cs:
            if c[0] == 'rat' and all(q.denominator == 1 for q in c[1]) and tuple(int(q) for q in c[1]) in S: continue
            L = sum(1 for pl in plist if contains(pl, c))
            if L > best[0]: best = (L, c)
    sizes = sorted((len(v) for v in planes.values()), reverse=True)
    return R, best, plist, planes, sizes

rows = [json.loads(l) for l in open('/home/user/Solving-hard-math-problems-/experiments/506/search-inversion/runs/results2.jsonl')]
seen = set(); summary = Counter()
for r in rows:
    u = r.get('universe', '')
    if not u.startswith('S2lat') or r.get('count') != r.get('formula') or r['n'] < 10: continue
    N = int(u.split('-')[1]); n = r['n']
    S = [tuple(int(t) for t in re.findall(r'-?\d+', lab)) for lab in r['S_labels']]
    key = (N, tuple(sorted(S)))
    if key in seen: continue
    seen.add(key)
    assert all(dot(p, p) == N for p in S) and len(set(S)) == n
    R, (L, c), plist, planes, sizes = best_centre(S, N)
    circ = R - L
    deg = sizes[0] == n
    note = ''
    if c is not None and c[0] == 'rat':
        img = [invert(p, c[1]) for p in S]
        nc, nl, csz, lsz = count_coplanar(img)
        note = f"projection count: circles={nc} lines={nl} maxcircle={csz[0]} lines={lsz} {'OK' if nc == circ else 'MISMATCH'}"
        antip = (csz[0] == n - 1 and csz[1:] == [3] * (nc - 1) and lsz == [3] * ((n - 1) // 2))
    else:
        note = "irrational centre (plane count only)"; antip = None
    ok = (circ == f(n)) and not deg
    summary[(ok, antip)] += 1
    print(f"{u} n={n}: claimed {r['count']} f={f(n)} | exact: R={R} L={L} circles={circ} degenerate={deg} centre={c if c is None or c[0]!='rat' else tuple(str(q) for q in c[1])} | {note} | antipodal-type={antip}", flush=True)
print("summary (claim reproduced, antipodal type):", dict(summary))
