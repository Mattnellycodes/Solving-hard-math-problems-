"""skel5.py n count_target : plan B for the single-5-block skeleton A = {0..4} (no other block of size >= 5).
At a point p of A the 4-blocks through p meet A in <= 1 further point; the derived structure at p is the 4-line
A - p plus d4(p) triples (pairwise sharing <= 1 point, each meeting A - p in <= 1 point).  Let d = max_{p in A} d4(p);
since at most one 4-block avoids A, b4 - 1 <= sum_{p in A} d4(p) <= 5 d.  Relabel so that d4(0) = d and enumerate the
derived structures at 0 up to isomorphism (canon with the 4-line as a distinguished block); CP-SAT completes each with
no further 4-block through 0 and d4(p) <= d on A."""
import sys, time, json, itertools
from math import comb
from pipe import *
n = int(sys.argv[1]); ctarget = int(sys.argv[2])
caps = caps_for(n, 'sg'); capD, capL = caps['capD'], caps['capL']
target = comb(n, 3) - ctarget
ellmax = capL // 3
b4min = -(-(target - 9 - ellmax) // 3)
dmax = (capD - 6) // 3
dmin = -(-(b4min - 1) // 5)
print(f'n={n} count<={ctarget}: need D+ell>={target}; single 5-block: b4>={b4min}, d in [{dmin},{dmax}]', flush=True)
A = frozenset(range(5)); m = n - 1
# derived structures at point 0: on points {1..n-1} (as 0..m-1 after shift by -1): 4-line {0,1,2,3} + triples
fourline = frozenset(range(4))
triples = [frozenset(T) for T in itertools.combinations(range(m), 3) if len(frozenset(T) & fourline) <= 1]
level = {None: [fourline]}; by_size = {}
k = 0
t0 = time.time()
while level and k < dmax:
    nxt = {}
    for code, fam in level.items():
        for T in triples:
            if T in fam or any(len(T & U) > 1 for U in fam[1:]): continue
            fam2 = fam + [T]
            c = canon(m, fam2)[0]
            if c not in nxt: nxt[c] = fam2
    k += 1; by_size[k] = list(nxt.values()); level = nxt
    print(f'  derived structures with {k} triples: {len(nxt)} classes [{time.time()-t0:.1f}s]', flush=True)
classes = {}
for d in range(dmin, dmax + 1):
    for j, P in enumerate(by_size.get(d, [])):
        fixed = [A] + [frozenset({0} | {v + 1 for v in T}) for T in P[1:]]
        t1 = time.time()
        st, sols = complete(n, fixed, target, caps, extra={'no4through': [0], 'max4': (list(A), d)}, time_limit=2400)
        F4s = {F4 for F4, L in sols}
        print(f'd={d} derived #{j}: status {st}, {len(sols)} (F4,L) solutions, {len(F4s)} distinct F4 [{time.time()-t1:.1f}s]', flush=True)
        assert st in ('OPTIMAL', 'INFEASIBLE'), st
        for F4 in F4s:
            F = fixed + sorted(F4, key=sorted)
            # must have no block of size >= 5 other than A: automatic (only 4-blocks added)
            code = canon(n, F)[0]
            classes.setdefault(code, dict(blocks=F, nlab=0))['nlab'] += 1
print(f'\n{len(classes)} isomorphism classes (single 5-block) with some completion [{time.time()-t0:.1f}s]', flush=True)
results = []
for k, (code, c) in enumerate(sorted(classes.items())):
    F = c['blocks']; D = sum(comb(len(B), 3) - 1 for B in F)
    Ls = line_sets(n, F, capL, target - D); ellmx = max((len(L) for L in Ls), default=-1)
    deg = sorted(sum(1 for B in F if p in B) for p in range(n)); cnt = comb(n, 3) - D - ellmx
    print(f'class {k}: sizes {sorted((len(B) for B in F), reverse=True)[:3]} b={len(F)} deg {deg} D={D} ell_max={ellmx} count={cnt} #linesets={len(Ls)}', flush=True)
    results.append(dict(idx=k, degrees=deg, D=D, ell_max=ellmx, count=cnt, n_line_sets=len(Ls), blocks=[sorted(B) for B in F],
                        line_sets=[sorted(sorted(l) for l in L) for L in Ls]))
json.dump(dict(n=n, ctarget=ctarget, results=results), open(f'out_skel5_n{n}_t{ctarget}.json', 'w'), indent=1)
