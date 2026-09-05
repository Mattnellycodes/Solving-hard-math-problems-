"""Driver: concurrency-event search (three blocks through one point) over symmetric families.
Usage: python3 search_concur.py m_min m_max [--nmax 12] [--out results_concur.json] [--budget 60]"""
import sys, os, json, math, time, argparse, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import symcore as S
from search import Family, Results, enumerate_families, nice_values, evaluate_params
from concur import concurrency_roots

HERE = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('m_min', type=int); ap.add_argument('m_max', type=int)
    ap.add_argument('--nmin', type=int, default=5)
    ap.add_argument('--nmax', type=int, default=12)
    ap.add_argument('--out', default=os.path.join(HERE, 'results_concur.json'))
    ap.add_argument('--budget', type=float, default=60.0)
    ap.add_argument('--G', type=int, default=160)
    ap.add_argument('--types', default='ABG')
    ap.add_argument('--max_orbits', type=int, default=3)
    args = ap.parse_args()
    results = Results(args.out)
    fams = enumerate_families(args.m_min, args.m_max, args.nmin, args.nmax, args.max_orbits, args.types)
    print(f"{len(fams)} families", flush=True)
    log = lambda s: print(s, flush=True)
    for m, orbs in fams:
        fam = Family(m, list(orbs))
        t0 = time.time()
        p = fam.p
        fixed = next((i for i, k in enumerate(fam.kinds) if k == 'r'), None)
        free = [i for i in range(p) if i != fixed]
        if not free:
            continue
        lists = [nice_values(fam.kinds[i], fam.m) for i in free]
        prods = list(itertools.product(*[range(min(len(L), 4)) for L in lists]))
        rng = np.random.default_rng(7)
        if len(prods) > 12: prods = [prods[j] for j in rng.choice(len(prods), 12, replace=False)]
        nroots = 0; nev = 0; seen = set()
        for pr in prods:
            if time.time() - t0 > args.budget: break
            base = np.ones(p)
            for i, j in zip(free, pr): base[i] = lists[free.index(i)][j]
            for i in free:
                if time.time() - t0 > args.budget: break
                try:
                    roots = concurrency_roots(fam, base, i, G=args.G)
                except Exception as ex:
                    log(f"  error {fam.name()} {ex}"); continue
                nroots += len(roots)
                for t in roots:
                    params = base.copy(); params[i] = t
                    key = tuple(np.round(params, 8))
                    if key in seen: continue
                    seen.add(key)
                    r = evaluate_params(fam, params); nev += 1
                    if r is None: continue
                    if results.add(fam, params, r, f'concur p{i}'):
                        log(f"  NEW n={fam.N} best={r['best']} (formula {S.formula(fam.N)}) {fam.name()} params={np.round(params,6).tolist()} nb={r['nblocks']} deg={r['deg']} lines={r['lines']} concur")
        log(f"{fam.name()} n={fam.N}: {nroots} concurrency roots, {nev} evals, {time.time()-t0:.1f}s")
    print("SUMMARY")
    for n in sorted(results.best):
        e = results.best[n][0]
        print(f"n={n:2d} formula={S.formula(n):3d} best={e['best']:3d} {e['family']} params={np.round(e['params'],5).tolist()}")
