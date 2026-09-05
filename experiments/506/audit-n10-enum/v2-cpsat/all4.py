"""all4.py n count_target : the all-4-block case (no rich block of size >= 5), own method.
Argument: D = 3 b4, ell <= floor(capL/3) =: ellmax, so b4 >= ceil((target - ellmax)/3); degrees d_p <= floor(capD/3);
a maximum-degree point p0 has d0 >= ceil(4 b4_min / n). Relabel p0 = 0; the blocks through 0 minus 0 form a
partial Steiner triple system (PSTS) with d0 triples on the other n-1 points (triples pairwise share <= 1 point).
For every isomorphism class of such PSTS (enumerated here by BFS + canon), CP-SAT enumerates all completions with
no further 4-block through 0 and all degrees <= d0. Every all-4-block structure is thus found (>= once)."""
import sys, time, json, itertools
from math import comb
from pipe import *

n = int(sys.argv[1]); ctarget = int(sys.argv[2])
caps = caps_for(n, 'sg'); capD, capL = caps['capD'], caps['capL']
target = comb(n, 3) - ctarget
ellmax = capL // 3
b4min = -(-(target - ellmax) // 3)
dmax = capD // 3
d0min = -(-4 * b4min // n)
print(f'n={n} count<={ctarget}: need D+ell>={target}; ellmax={ellmax}, b4>={b4min}, degrees<={dmax}, max degree d0 in [{d0min},{dmax}]', flush=True)
m = n - 1
t0 = time.time()
# PSTS(m) classes by number of triples
triples = [frozenset(T) for T in itertools.combinations(range(m), 3)]
level = {(): []}; psts_by_size = {0: [[]]}
k = 0
while level and k < dmax:
    nxt = {}
    for code, fam in level.items():
        for T in triples:
            if T in fam or any(len(T & U) > 1 for U in fam): continue
            fam2 = fam + [T]
            c = canon(m, fam2)[0]
            if c not in nxt: nxt[c] = fam2
    k += 1
    psts_by_size[k] = list(nxt.values())
    level = nxt
    print(f'  PSTS({m}) with {k} triples: {len(nxt)} classes', flush=True)
classes = {}
for d0 in range(d0min, dmax + 1):
    for j, P in enumerate(psts_by_size.get(d0, [])):
        fixed = [frozenset({0} | {v + 1 for v in T}) for T in P]
        t1 = time.time()
        st, sols = complete(n, fixed, target, caps, extra={'no4through': [0], 'maxdeg': d0}, time_limit=2400)
        F4s = {F4 for F4, L in sols}
        print(f'd0={d0} PSTS #{j}: status {st}, {len(sols)} (F4,L) solutions, {len(F4s)} distinct F4 [{time.time()-t1:.1f}s]', flush=True)
        assert st in ('OPTIMAL', 'INFEASIBLE'), st
        for F4 in F4s:
            F = fixed + sorted(F4, key=sorted)
            code = canon(n, F)[0]
            classes.setdefault(code, dict(blocks=F, nlab=0))['nlab'] += 1
print(f'\n{len(classes)} isomorphism classes of all-4-block structures with some completion [{time.time()-t0:.1f}s]', flush=True)
results = []
for k, (code, c) in enumerate(sorted(classes.items())):
    F = c['blocks']
    D = sum(comb(len(B), 3) - 1 for B in F)
    Ls = line_sets(n, F, capL, target - D)
    ellmx = max(len(L) for L in Ls) if Ls else -1
    deg = sorted(sum(1 for B in F if p in B) for p in range(n))
    auto = aut_group_order(n, F)[0]
    cnt = comb(n, 3) - D - ellmx
    surv = {'A': 0, 'AB': 0, 'ABK': 0, 'ABKC': 0}
    for L in Ls:
        f = full_check(n, F, L)
        okA = not f['A']; okB = okA and not f['B']; okK = okB and not f['K']; okC = okK and not f['C']
        surv['A'] += okA; surv['AB'] += okB; surv['ABK'] += okK; surv['ABKC'] += okC
    print(f'class {k}: b4={len(F)} deg {deg} |Aut| {auto} D={D} ell_max={ellmx} count={cnt} #linesets={len(Ls)} survivors(A/AB/ABK/ABKC)={surv}', flush=True)
    print('   blocks', [sorted(B) for B in F])
    results.append(dict(idx=k, degrees=deg, aut=auto, D=D, ell_max=ellmx, count=cnt, n_line_sets=len(Ls), survivors=surv,
                        blocks=[sorted(B) for B in F], line_sets=[sorted(sorted(l) for l in L) for L in Ls]))
json.dump(dict(n=n, ctarget=ctarget, results=results), open(f'out_all4_n{n}_t{ctarget}' + ('_noQ' if not USE_Q else '') + '.json', 'w'), indent=1)
print(f'total time {time.time()-t0:.1f}s')
