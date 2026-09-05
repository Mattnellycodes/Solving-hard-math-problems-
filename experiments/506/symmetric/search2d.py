"""2-parameter search: trace the curve {condition A = 0} in the plane of two free parameters
(p_i, p_j) and locate the points on it where a second condition vanishes (simultaneous
coincidences).  Every such point is evaluated with the inversion-centre optimisation.

Usage: python3 search2d.py m_min m_max [--nmax 16] [--out results_2d.json] [--budget 150]
"""
import sys, os, json, math, time, argparse, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import symcore as S
from search import Family, Results, enumerate_families, nice_values, evaluate_params

HERE = os.path.dirname(os.path.abspath(__file__))


def roots_1d_all(fam, base, i, G):
    """dict: condition index -> sorted list of roots in parameter i (others at base)."""
    base = np.asarray(base, float)
    grid = fam.grid(fam.kinds[i], G)
    pb = np.repeat(base[None, :], G, 0); pb[:, i] = grid
    V = fam.conds_batch(pb)
    sg = np.sign(V)
    ident = np.abs(V).max(0) < 1e-9
    ch = np.where((sg[:-1] * sg[1:]) < 0)
    lo = grid[ch[0]].copy(); hi = grid[ch[0] + 1].copy(); cols = ch[1]
    keep = ~ident[cols]
    lo, hi, cols = lo[keep], hi[keep], cols[keep]
    if len(cols) == 0: return {}, V, grid
    flo = V[ch[0][keep], cols].copy()
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        pb = np.repeat(base[None, :], len(mid), 0); pb[:, i] = mid
        fm = fam.cond_single(pb, cols)
        left = (np.sign(fm) == np.sign(flo)) | (fm == 0)
        lo = np.where(left, mid, lo); flo = np.where(left, fm, flo); hi = np.where(left, hi, mid)
        if np.all(hi - lo < 1e-11): break
    out = {}
    for c, t in zip(cols, 0.5 * (lo + hi)):
        out.setdefault(int(c), []).append(float(t))
    for c in out: out[c].sort()
    return out, V, grid


def solve_on_curve(fam, base, i, j, c, pj, ti_guess, width):
    """solve condition c = 0 for p_i near ti_guess at p_j = pj. Returns t or None."""
    base = np.asarray(base, float).copy(); base[j] = pj
    lo, hi = ti_guess - width, ti_guess + width
    if fam.kinds[i] == 'a':
        lo = max(lo, 1e-4); hi = min(hi, math.pi / fam.m - 1e-4)
    else:
        lo = max(lo, 1e-3)
    pb = np.repeat(base[None, :], 2, 0); pb[0, i] = lo; pb[1, i] = hi
    f = fam.cond_single(pb, [c, c])
    if np.sign(f[0]) == np.sign(f[1]):
        # try a finer scan inside
        ts = np.linspace(lo, hi, 25)
        pb = np.repeat(base[None, :], 25, 0); pb[:, i] = ts
        v = fam.cond_single(pb, [c] * 25)
        s = np.sign(v)
        k = np.where(s[:-1] * s[1:] < 0)[0]
        if len(k) == 0: return None
        # nearest to guess
        k = k[np.argmin(np.abs(ts[k] - ti_guess))]
        lo, hi = ts[k], ts[k + 1]; f = np.array([v[k], v[k + 1]])
    flo = f[0]
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        pb = base[None, :].copy(); pb[0, i] = mid
        fm = fam.cond_single(pb, [c])[0]
        if np.sign(fm) == np.sign(flo) or fm == 0: lo, flo = mid, fm
        else: hi = mid
        if hi - lo < 1e-12: break
    return 0.5 * (lo + hi)


def trace_pairs(fam, results, base, i, j, G1=90, G2=250, time_budget=150.0, log=print):
    """For parameter pair (i, j): trace all condition curves over a grid in p_j and find crossings."""
    t0 = time.time()
    gridj = fam.grid(fam.kinds[j], G1)
    data = []       # per grid point: (roots dict, V (G2,K), grid_i)
    for pj in gridj:
        if time.time() - t0 > time_budget: break
        b = np.asarray(base, float).copy(); b[j] = pj
        data.append(roots_1d_all(fam, b, i, G2))
    found = 0; evals = 0
    K = fam.K
    cand = []
    for g in range(len(data) - 1):
        R0, _, _ = data[g]; R1, _, _ = data[g + 1]
        for c in set(R0) & set(R1):
            for t0_ in R0[c]:
                # match nearest root at next grid point
                t1 = min(R1[c], key=lambda t: abs(t - t0_))
                step = max(abs(t1 - t0_), 1e-6)
                if abs(t1 - t0_) > 0.3 * max(abs(t0_), 0.1) + 0.05: continue
                # evaluate all conditions at both curve points
                b0 = np.asarray(base, float).copy(); b0[j] = gridj[g]; b0[i] = t0_
                b1 = np.asarray(base, float).copy(); b1[j] = gridj[g + 1]; b1[i] = t1
                V = fam.conds_batch(np.stack([b0, b1]))
                s = np.sign(V)
                cross = np.where((s[0] * s[1] < 0) & (np.abs(V).max(0) > 1e-9))[0]
                for c2 in cross:
                    if c2 == c: continue
                    cand.append((c, int(c2), gridj[g], gridj[g + 1], t0_, t1, V[0, c2], V[1, c2]))
    log(f"   pair ({i},{j}): {len(data)} grid pts, {len(cand)} crossing candidates ({time.time()-t0:.0f}s)")
    seenpts = set()
    def F(pi_, pj_, c, c2):
        b = np.asarray(base, float).copy(); b[i] = pi_; b[j] = pj_
        return fam.cond_single(np.repeat(b[None, :], 2, 0), [c, c2])
    for (c, c2, pj0, pj1, ti0, ti1, f0, f1) in cand:
        if time.time() - t0 > 2 * time_budget: break
        # 2-D Newton with finite-difference Jacobian from the segment midpoint
        x = np.array([0.5 * (ti0 + ti1), 0.5 * (pj0 + pj1)])
        box = (min(ti0, ti1) - 0.1 * abs(ti1 - ti0) - 1e-6, max(ti0, ti1) + 0.1 * abs(ti1 - ti0) + 1e-6, pj0 - 1e-9, pj1 + 1e-9)
        ok = False
        for it in range(25):
            f = F(x[0], x[1], c, c2)
            if np.abs(f).max() < 1e-13: ok = True; break
            h = 1e-7
            J = np.empty((2, 2))
            J[:, 0] = (F(x[0] + h, x[1], c, c2) - f) / h
            J[:, 1] = (F(x[0], x[1] + h, c, c2) - f) / h
            try:
                dx = np.linalg.solve(J, -f)
            except np.linalg.LinAlgError:
                break
            x = x + dx
            if not (box[0] - 0.05 <= x[0] <= box[1] + 0.05 and box[2] - 0.05 * (pj1 - pj0) - 1e-3 <= x[1] <= box[3] + 0.05 * (pj1 - pj0) + 1e-3): break
        if not ok:
            f = F(x[0], x[1], c, c2)
            if np.abs(f).max() > 1e-10: continue
        params = np.asarray(base, float).copy(); params[i] = x[0]; params[j] = x[1]
        key = tuple(np.round(params, 7))
        if key in seenpts: continue
        seenpts.add(key)
        r = evaluate_params(fam, params); evals += 1
        if r is None: continue
        if results.add(fam, params, r, f'2d p{i}p{j} conds {c},{c2}'):
            log(f"  NEW n={fam.N} best={r['best']} (formula {S.formula(fam.N)}) {fam.name()} params={np.round(params,6).tolist()} nb={r['nblocks']} deg={r['deg']} lines={r['lines']}")
            found += 1
    log(f"   pair ({i},{j}): {evals} evaluations, {found} new, {time.time()-t0:.0f}s")


def search_family_2d(fam, results, budget, log=print):
    p = fam.p
    fixed = next((i for i, k in enumerate(fam.kinds) if k == 'r'), None)
    free = [i for i in range(p) if i != fixed]
    if len(free) < 2: return
    base = np.ones(p)
    for i in free:
        base[i] = nice_values(fam.kinds[i], fam.m)[0]
    others = [i for i in free]
    for i, j in itertools.combinations(free, 2):
        rest = [k for k in free if k not in (i, j)]
        # other free params at a couple of nice values
        combos = [()] if not rest else list(itertools.product(*[nice_values(fam.kinds[k], fam.m)[:2] for k in rest]))
        for combo in combos[:2]:
            b = base.copy()
            for k, v in zip(rest, combo): b[k] = v
            trace_pairs(fam, results, b, i, j, time_budget=budget, log=log)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('m_min', type=int); ap.add_argument('m_max', type=int)
    ap.add_argument('--nmin', type=int, default=5)
    ap.add_argument('--nmax', type=int, default=16)
    ap.add_argument('--out', default=os.path.join(HERE, 'results_2d.json'))
    ap.add_argument('--budget', type=float, default=150.0)
    ap.add_argument('--max_orbits', type=int, default=3)
    args = ap.parse_args()
    results = Results(args.out)
    fams = enumerate_families(args.m_min, args.m_max, args.nmin, args.nmax, args.max_orbits)
    fams = [(m, o) for m, o in fams if S.nparams(list(o)) >= 3]
    print(f"{len(fams)} families with >= 2 free parameters", flush=True)
    for m, orbs in fams:
        fam = Family(m, list(orbs))
        t0 = time.time()
        print(f"{fam.name()} n={fam.N} p={fam.p} K={fam.K}", flush=True)
        search_family_2d(fam, results, args.budget, log=lambda s: print(s, flush=True))
    print("SUMMARY")
    for n in sorted(results.best):
        e = results.best[n][0]
        print(f"n={n:2d} formula={S.formula(n):3d} best={e['best']:3d} {e['family']} params={np.round(e['params'],5).tolist()}")
