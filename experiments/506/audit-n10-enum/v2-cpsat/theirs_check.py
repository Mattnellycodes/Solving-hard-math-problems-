"""Independent re-check of the 9 structures listed in ../../n10-enum/orderly3/n10_sg_t32.json (data only):
own line-set enumeration (with and without the Q_p cap), own hereditary checks, |Aut|, canonical codes."""
import json, sys, os
from math import comb
import pipe
from pipe import *
their = json.load(open('/home/user/Solving-hard-math-problems-/experiments/506/n10-enum/orderly3/n10_sg_t32.json'))
n = 10; caps = caps_for(n); target = 88
out = []
for i, r in enumerate(their['results']):
    F = [frozenset(B) for B in r['blocks']]
    D = sum(comb(len(B), 3) - 1 for B in F)
    code = canon(n, F)[0]
    auto = aut_group_order(n, F)[0]
    pipe.USE_Q = False
    Ls_noQ = line_sets(n, F, caps['capL'], target - D)
    pipe.USE_Q = True
    Ls = line_sets(n, F, caps['capL'], target - D)
    surv = {'A': [], 'AB': [], 'ABK': [], 'ABKC': []}
    for L in Ls_noQ:
        f = full_check(n, F, L)
        okA = not f['A']; okB = okA and not f['B']; okK = okB and not f['K']; okC = okK and not f['C']
        for key, ok in (('A', okA), ('AB', okB), ('ABK', okK), ('ABKC', okC)):
            if ok: surv[key].append(sorted(sorted(l) for l in L))
    ellmax = max((len(L) for L in Ls_noQ), default=-1)
    ellmaxQ = max((len(L) for L in Ls), default=-1)
    print(f"their #{i}: sizes {sorted((len(B) for B in F), reverse=True)[:5]}... b={len(F)} deg {sorted(sum(1 for B in F if p in B) for p in range(n))} |Aut|={auto} D={D} "
          f"ell_max={ellmax} (their {r['ell_max']}) count={comb(n,3)-D-ellmax} (their {r['count_min']}) #linesets noQ={len(Ls_noQ)} (their {len(r['line_sets'])}) withQ={len(Ls)} "
          f"survivors A/AB/ABK/ABKC = {[len(surv[k]) for k in ('A','AB','ABK','ABKC')]}", flush=True)
    for k in ('ABKC', 'ABK'):
        if surv[k] and k == 'ABKC':
            for L in surv[k]: print('     survives all tiers, count', comb(n,3)-D-len(L), 'lines', L)
    out.append(dict(idx=i, code=list(code), aut=auto, D=D, ell_max=ellmax, n_line_sets_noQ=len(Ls_noQ), n_line_sets_Q=len(Ls),
                    survivors={k: len(v) for k, v in surv.items()}, surviving_ABK=surv['ABK'], surviving_ABKC=surv['ABKC']))
json.dump(out, open('theirs_check.json', 'w'), indent=1)
