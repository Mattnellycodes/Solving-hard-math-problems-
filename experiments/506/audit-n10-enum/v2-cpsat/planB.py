"""planB.py n count_target skel_index [--orbit-point q]: generic decomposition for one skeleton class S (same
indexing as run.py: skeleton classes from enumerate_skeletons, nonempty ones in order).
Choose a point q (default: a point of the first big block) and let O = orbit of q under Aut(S).  Let
d = max_{p in O} d4(p).  Relabelling by an element of Aut(S) we may assume d4(q) = d.  The 4-blocks through q form a
'derived structure' (triples on P - q, pairwise sharing <= 1 point, each meeting B - q in <= 1 point for big B
through q and B in <= 2 points for big B not through q).  These are enumerated up to Aut(S)_q-equivalence (canonical
form of S + the 4-blocks through q with q individualised).  CP-SAT then completes each with no further 4-block
through q and d4 <= d on O.  Every completion of S is found at least once."""
import sys, time, json, itertools
from math import comb
from pipe import *
from canon import canon, aut_group_order
n = int(sys.argv[1]); ctarget = int(sys.argv[2]); si = int(sys.argv[3])
caps = caps_for(n, 'sg'); capD, capL = caps['capD'], caps['capL']
target = comb(n, 3) - ctarget
skels = [s for s in enumerate_skeletons(n, capD, verbose=False) if s]
S = skels[si]
order, G = aut_group_order(n, S)
q = int(sys.argv[sys.argv.index('--orbit-point') + 1]) if '--orbit-point' in sys.argv else min(sorted(S, key=lambda B: (len(B), sorted(B)))[0])
O = sorted({g[q] for g in G})
print(f'n={n} count<={ctarget} need D+ell>={target}; skeleton {si} = {[sorted(B) for B in S]} |Aut|={order}; q={q}, orbit O={O}', flush=True)
used_q = sum(comb(len(B) - 1, 2) for B in S if q in B)
dmax = (capD - used_q) // 3
others = [p for p in range(n) if p != q]
cand = []
for T in itertools.combinations(others, 3):
    T = frozenset(T)
    if all((len(T & B) <= 1) if q in B else (len(T & B) <= 2) for B in S): cand.append(T)
init = [[q], others]
level = {None: []}; by_size = {0: [[]]}
t0 = time.time(); k = 0
while level and k < dmax:
    nxt = {}
    for code, fam in level.items():
        for T in cand:
            if T in fam or any(len(T & U) > 1 for U in fam): continue
            fam2 = fam + [T]
            c = canon(n, list(S) + [T2 | {q} for T2 in fam2], init_cells=init)[0]
            if c not in nxt: nxt[c] = fam2
    k += 1; by_size[k] = list(nxt.values()); level = nxt
    print(f'  derived structures at q with {k} triples: {len(nxt)} classes [{time.time()-t0:.1f}s]', flush=True)
classes = {}
for d in range(0, dmax + 1):
    for j, P in enumerate(by_size.get(d, [])):
        fixed = list(S) + [T | {q} for T in P]
        t1 = time.time()
        st, sols = complete(n, fixed, target, caps, extra={'no4through': [q], 'max4': (O, d)}, time_limit=2400)
        F4s = {F4 for F4, L in sols}
        print(f'd={d} derived #{j}: status {st}, {len(sols)} (F4,L) solutions, {len(F4s)} distinct F4 [{time.time()-t1:.1f}s]', flush=True)
        assert st in ('OPTIMAL', 'INFEASIBLE'), st
        for F4 in F4s:
            F = fixed + sorted(F4, key=sorted)
            code = canon(n, F)[0]
            classes.setdefault(code, dict(blocks=F, nlab=0))['nlab'] += 1
print(f'\n{len(classes)} isomorphism classes (skeleton {si}) with some completion [{time.time()-t0:.1f}s]', flush=True)
results = []
for k, (code, c) in enumerate(sorted(classes.items())):
    F = c['blocks']; D = sum(comb(len(B), 3) - 1 for B in F)
    Ls = line_sets(n, F, capL, target - D); ellmx = max((len(L) for L in Ls), default=-1)
    deg = sorted(sum(1 for B in F if p in B) for p in range(n)); cnt = comb(n, 3) - D - ellmx
    auto = aut_group_order(n, F)[0]
    surv = {'A': 0, 'AB': 0, 'ABK': 0, 'ABKC': 0}
    for L in Ls:
        f = full_check(n, F, L)
        okA = not f['A']; okB = okA and not f['B']; okK = okB and not f['K']; okC = okK and not f['C']
        surv['A'] += okA; surv['AB'] += okB; surv['ABK'] += okK; surv['ABKC'] += okC
    print(f'class {k}: sizes {sorted((len(B) for B in F), reverse=True)[:4]} b={len(F)} deg {deg} |Aut| {auto} D={D} ell_max={ellmx} count={cnt} #linesets={len(Ls)} survivors={surv}', flush=True)
    print('   blocks', [sorted(B) for B in F])
    results.append(dict(idx=k, skeleton=si, degrees=deg, aut=auto, D=D, ell_max=ellmx, count=cnt, n_line_sets=len(Ls), survivors=surv,
                        blocks=[sorted(B) for B in F], line_sets=[sorted(sorted(l) for l in L) for L in Ls]))
json.dump(dict(n=n, ctarget=ctarget, skeleton=si, results=results), open(f'out_planB_n{n}_t{ctarget}_s{si}' + ('_noQ' if not USE_Q else '') + '.json', 'w'), indent=1)
