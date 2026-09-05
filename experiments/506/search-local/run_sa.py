"""Simulated-annealing / tabu subset search over rich universes (Erdős #506, search-local agent).

usage: python3 run_sa.py <universe names or 'all'> [--ns 9-16] [--iters 300000] [--restarts 6] [--seed 0]
Appends one JSON line per (universe, n) to runs/sa_results.jsonl; records with true count <= f(n)
are appended (with exact coordinates when available) to runs/sa_records.jsonl.
"""
import sys, time, json, argparse, math
import numpy as np
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/search-local')
from engine import *
import universes as UV
import sympy as sp
from fractions import Fraction as Fr

HERE = '/home/user/Solving-hard-math-problems-/experiments/506/search-local'


def catalogue():
    C = {}
    C['grid4'] = lambda: UV.grid(4)
    C['grid5'] = lambda: UV.grid(5)
    C['grid6'] = lambda: UV.grid(6)
    C['grid7'] = lambda: UV.grid(7)
    C['rect4x6'] = lambda: UV.rect(4, 6)
    C['rect3x7'] = lambda: UV.rect(3, 7)
    C['grid5inv-c'] = lambda: UV.grid_inversions(5, [(2, 2)], [1, 2, 4, 5, 8], name='grid5inv-c')
    C['grid5inv-corner'] = lambda: UV.grid_inversions(5, [(0, 0)], [1, 2, 4, 5, 8, 10], name='grid5inv-corner')
    C['grid4inv-multi'] = lambda: UV.grid_inversions(4, [(0, 0), (1, 1), (3, 3), (1, 0)], [1, 2, 4, 5], name='grid4inv-multi')
    C['grid5inv-multi'] = lambda: UV.grid_inversions(5, [(2, 2), (0, 0), (1, 1), (2, 0)], [1, 2, 4, 5], name='grid5inv-multi')
    C['grid6inv-c'] = lambda: UV.grid_inversions(6, [(2, 2), (3, 3)], [1, 2, 4, 5], name='grid6inv-c')
    C['trilat3'] = lambda: UV.tri_lattice(3)
    C['trilat4'] = lambda: UV.tri_lattice(4)
    C['trilat5'] = lambda: UV.tri_lattice(5)
    phi = (1 + sp.sqrt(5)) / 2
    s2, s3 = sp.sqrt(2), sp.sqrt(3)
    C['poly5-phi'] = lambda: UV.polygons(5, [1, phi, phi**2, phi**3, 1/phi, 1/phi**2, 1, phi, phi**2, 1/phi], [0]*6 + [1]*4, name='poly5-phi')
    C['poly8-s2'] = lambda: UV.polygons(8, [1, s2, 2, 2*s2, 4, 1+s2, 1, s2, 2, 1+s2, s2-1], [0]*6 + [1]*5, name='poly8-s2')
    C['poly6-s3'] = lambda: UV.polygons(6, [1, s3, 2, 3, 2*s3, sp.sqrt(7), 1, s3, 2, 3, sp.sqrt(7)], [0]*6 + [1]*5, name='poly6-s3')
    C['poly12'] = lambda: UV.polygons(12, [1, 2, s3, 2+s3, 2-s3, (sp.sqrt(6)+s2)/2, 1, 2, s3, 2+s3, 2-s3], [0]*6 + [1]*5, name='poly12')
    C['poly7'] = lambda: UV.polygons(7, [1, 2, 1/(2*sp.cos(sp.pi/7)), 2*sp.cos(sp.pi/7), 2*sp.cos(2*sp.pi/7), 1, 2, 1/(2*sp.cos(sp.pi/7)), 2*sp.cos(sp.pi/7), 2*sp.cos(2*sp.pi/7)], [0]*5 + [1]*5, name='poly7')
    C['poly9'] = lambda: UV.polygons(9, [1, 2, 2*sp.cos(sp.pi/9), 2*sp.cos(2*sp.pi/9), 1/(2*sp.cos(sp.pi/9)), 1, 2, 2*sp.cos(sp.pi/9), 2*sp.cos(2*sp.pi/9)], [0]*5 + [1]*4, name='poly9')
    C['poly10'] = lambda: UV.polygons(10, [1, phi, phi**2, 1/phi, 2, 1, phi, phi**2, 1/phi, 2], [0]*5 + [1]*5, name='poly10')
    C['poly14'] = lambda: UV.polygons(14, [1, 2, 2*sp.cos(sp.pi/7), 2*sp.cos(2*sp.pi/7), 1, 2, 2*sp.cos(sp.pi/7)], [0]*4 + [1]*3, name='poly14')
    C['poly16'] = lambda: UV.polygons(16, [1, s2, 2, 1+s2, 1, s2, 2], [0]*4 + [1]*3, name='poly16')
    C['poly18'] = lambda: UV.polygons(18, [1, 2, 2*sp.cos(sp.pi/9), 1, 2, 2*sp.cos(sp.pi/9)], [0]*3 + [1]*3, name='poly18')
    C['poly20'] = lambda: UV.polygons(20, [1, phi, 2, 1, phi, 2], [0]*3 + [1]*3, name='poly20')
    C['poly24'] = lambda: UV.polygons(24, [1, 2, s3, 1, 2, s3], [0]*3 + [1]*3, name='poly24')
    C['poly30'] = lambda: UV.polygons(30, [1, phi, 2, 1, phi], [0]*3 + [1]*2, name='poly30')
    # orthocentric closures (rational)
    C['ortho-A1'] = lambda: UV.ortho_closure([(0, 3), (1, 0), (3, 0)], levels=1, name='ortho-A1')
    C['ortho-A2'] = lambda: UV.ortho_closure([(0, 3), (1, 0), (3, 0)], levels=2, cap=250, name='ortho-A2')
    C['ortho-B2'] = lambda: UV.ortho_closure([(0, 0), (4, 0), (1, 3)], levels=2, cap=250, name='ortho-B2')
    C['ortho-C2'] = lambda: UV.ortho_closure([(0, 0), (6, 0), (2, 4)], levels=2, cap=250, name='ortho-C2')
    C['ortho-D2'] = lambda: UV.ortho_closure([(0, 0), (20, 0), (5, 15)], levels=2, cap=250, name='ortho-D2')
    C['ortho-E2'] = lambda: UV.ortho_closure([(0, 0), (12, 0), (3, 9), (6, 3)], levels=2, cap=250, name='ortho-E2')
    C['ortho-F1'] = lambda: UV.ortho_closure([(0, 0), (4, 0), (0, 4), (4, 4), (2, 1)], levels=1, cap=300, name='ortho-F1')
    C['ortho-G2'] = lambda: UV.ortho_closure([(0, 0), (5, 0), (5, 5), (0, 5)], levels=2, cap=250, name='ortho-G2')
    C['ortho-H2'] = lambda: UV.ortho_closure([(0, 0), (8, 0), (2, 4)], levels=2, cap=250, ops=('H', 'O', 'M', 'F', 'RH', 'RM'), name='ortho-H2')
    # sphere universes
    C['icosa-I-v'] = lambda: UV.icosa_family(('I',), 'vertex')
    C['icosa-ID-v'] = lambda: UV.icosa_family(('I', 'D'), 'vertex')
    C['icosa-ID-f'] = lambda: UV.icosa_family(('I', 'D'), 'face')
    C['icosa-ID-g'] = lambda: UV.icosa_family(('I', 'D'), 'generic')
    C['icosa-IDE-v'] = lambda: UV.icosa_family(('I', 'D', 'E'), 'vertex')
    C['icosa-IDE-g'] = lambda: UV.icosa_family(('I', 'D', 'E'), 'generic')
    C['octa-123-v'] = lambda: UV.octa_shells((1, 2, 3), 'vertex')
    C['octa-123-g'] = lambda: UV.octa_shells((1, 2, 3), 'generic')
    C['octa-12356-v'] = lambda: UV.octa_shells((1, 2, 3, 5, 6), 'vertex')
    C['octa-12356-e'] = lambda: UV.octa_shells((1, 2, 3, 5, 6), 'edge')
    C['octa-123569-v'] = lambda: UV.octa_shells((1, 2, 3, 5, 6, 9), 'vertex')
    # hyperbola "circle orchard": xy = 1 at t = ±2^k, plus centre, infinity
    def hyp():
        pts, lab = [], []
        for k in range(-3, 4):
            for s in (1, -1):
                t = Fr(s) * Fr(2) ** k
                pts.append((t, 1 / t)); lab.append(f"hyp(t={t})")
        pts.append((Fr(0), Fr(0))); lab.append('C')
        return UV._from_exact('hyperbola14', pts, lab, True)
    C['hyperbola14'] = hyp

    def hyp3():
        pts, lab = [], []
        for k in range(-2, 3):
            for s in (1, -1):
                for base in (2, 3):
                    t = Fr(s) * Fr(base) ** k
                    if (t, 1 / t) not in pts:
                        pts.append((t, 1 / t)); lab.append(f"hyp(t={t})")
        pts.append((Fr(0), Fr(0))); lab.append('C')
        return UV._from_exact('hyperbola23', pts, lab, True)
    C['hyperbola23'] = hyp3
    # conics
    C['ellipse24-2-1'] = lambda: UV.ellipse(24, 2, 1, ('C',))
    C['ellipse24-2-1F'] = lambda: UV.ellipse(24, 2, 1, ('C', 'F'))
    C['ellipse30-3-2'] = lambda: UV.ellipse(30, 3, 2, ('C',))
    C['ellipse20-5-3F'] = lambda: UV.ellipse(20, 5, 3, ('C', 'F'))
    C['ellipse24-2-1-rot'] = lambda: UV.ellipse(24, 2, 1, ('C',), theta0=sp.pi / 24, name='ellipse24-2-1-rot')
    C['hyperbola-pow'] = lambda: UV.hyperbola_exact([Fr(s) * Fr(b) ** k for b in (2, 3) for k in range(-2, 3) for s in (1, -1)] + [Fr(3, 2), Fr(2, 3), Fr(-3, 2), Fr(-2, 3)], name='hyperbola-pow')
    # polygon diagonal closures (selectable concurrency points)
    C['pdiag8'] = lambda: UV.polygon_diag_closure(8, min_mult=3, max_pts=120)
    C['pdiag10'] = lambda: UV.polygon_diag_closure(10, min_mult=3, max_pts=120)
    C['pdiag12'] = lambda: UV.polygon_diag_closure(12, min_mult=3, max_pts=120)
    C['pdiag16'] = lambda: UV.polygon_diag_closure(16, min_mult=4, max_pts=120)
    C['pdiag18'] = lambda: UV.polygon_diag_closure(18, min_mult=4, max_pts=120)
    C['pdiag24'] = lambda: UV.polygon_diag_closure(24, min_mult=5, max_pts=120)
    C['pdiag30'] = lambda: UV.polygon_diag_closure(30, min_mult=5, max_pts=120)
    C['pdiag8x2'] = lambda: UV.polygon_diag_closure(8, min_mult=4, max_pts=100, radius2=1 + sp.sqrt(2))
    C['pdiag12x2'] = lambda: UV.polygon_diag_closure(12, min_mult=5, max_pts=100, radius2=2 + sp.sqrt(3))
    return C


def closure_universe(U, name, min_mult=3, max_pts=150):
    """add the highest-multiplicity arrangement points of U as selectable points."""
    ext = UV.multiplicity_points(U, min_mult=min_mult, max_pairs=4_000_000)
    ext.sort(key=lambda e: -e[2])
    ext = ext[:max_pts]
    P = np.vstack([U.P, np.array([[e[0], e[1]] for e in ext])]) if ext else U.P
    labels = list(U.labels or [U.label(i) for i in range(U.N)]) + [f"X{e[2]}#{i}" for i, e in enumerate(ext)]
    return Universe(name, P, has_inf=U.has_inf, exact=None, labels=labels)


def exact_coords(U, S):
    """exact coordinate strings for S (None if unavailable)"""
    if U.exact is None and hasattr(U, 'exact_lazy'):
        out = []
        for i in S:
            i = int(i)
            if U.has_inf and i == U.N:
                out.append('inf')
            else:
                x, y = UV.resolve_lazy_exact(U, i)
                out.append([str(x), str(y)])
        return out
    if U.exact is None:
        return None
    out = []
    for i in S:
        i = int(i)
        if U.has_inf and i == U.N:
            out.append('inf')
        else:
            x, y = U.exact[i]
            out.append([str(x), str(y)])
    return out


def done_keys():
    keys = set()
    try:
        for line in open(f"{HERE}/runs/sa_results.jsonl"):
            try:
                r = json.loads(line); keys.add((r['universe'], r['n']))
            except Exception:
                pass
    except FileNotFoundError:
        pass
    return keys


def run_universe(U, ns, iters, restarts, seed, log, t_limit=None, tabu=4, bias=4, skip_done=False, budget=4e8):
    res = {}
    rng = np.random.default_rng(seed)
    avg_deg = U.pptr[-1] / max(U.M, 1)          # average number of blocks through a point
    use_table = U.M <= 320
    if use_table:
        tab = build_triple_table(U)
        iters_eff = iters
        log(f"  avg blocks/point={avg_deg:.0f}; triple-table kernel; iters per restart={iters_eff}")
    else:
        iters_eff = int(min(iters, max(50_000, budget / avg_deg)))
        log(f"  avg blocks/point={avg_deg:.0f}; CSR kernel; iters per restart={iters_eff}")
    done = done_keys() if skip_done else set()
    prev_best = None
    for n in ns:
        if n >= U.M:
            continue
        if (U.name, n) in done:
            prev_best = None
            continue
        iters = iters_eff
        t = time.time()
        best = (BIG, None, None)
        allbest = []
        for r in range(restarts):
            # random non-degenerate init; restart 0 is seeded from the best (n-1)-subset plus a point
            init = None
            if r == 0 and prev_best is not None and len(prev_best) == n - 1:
                for _ in range(100):
                    extra = int(rng.integers(U.M))
                    if extra in set(int(i) for i in prev_best):
                        continue
                    cand = np.array(list(prev_best) + [extra], dtype=np.int64)
                    v, o, nb = k_eval(cand, n, U.T, U.bptr, U.bmem, U.pptr, U.pblk)
                    if v < BIG:
                        init = cand; break
            if init is None:
                for _ in range(100):
                    init = rng.choice(U.M, size=n, replace=False).astype(np.int64)
                    v, o, nb = k_eval(init, n, U.T, U.bptr, U.bmem, U.pptr, U.pblk)
                    if v < BIG:
                        break
            T0 = 3.0 if r % 2 == 0 else 1.5
            if use_table:
                b, S, O = t_sa(n, U.M, U.T, tab, U.bptr, U.bmem, len(U.blocks), iters, T0, 0.15, int(rng.integers(1 << 30)), init, tabu, bias)
            else:
                b, S, O = k_sa(n, U.M, U.T, U.bptr, U.bmem, U.pptr, U.pblk, iters, T0, 0.15, int(rng.integers(1 << 30)), init, tabu, bias)
            allbest.append((int(b), S.copy(), int(O)))
            if b < best[0]:
                best = (int(b), S.copy(), int(O))
            if t_limit and time.time() - t > t_limit:
                break
        # true count on the distinct best few subsets
        allbest.sort(key=lambda z: z[0])
        seen = set(); evald = []
        for b, S, O in allbest[:5]:
            key = tuple(sorted(int(i) for i in S))
            if key in seen:
                continue
            seen.add(key)
            tc, Ot, d, nb, lines = true_count(U, S)
            evald.append((tc, b, S, Ot, d, nb, lines))
        evald.sort(key=lambda z: z[0])
        tc, b, S, Ot, d, nb, lines = evald[0]
        if tc < b:
            # the proxy is exact; a smaller float value needs an exact certificate
            ok = False
            try:
                from verify_exact import certify
                ex = exact_coords(U, S)
                if ex is not None:
                    import signal
                    def _h(signum, frame):
                        raise TimeoutError()
                    signal.signal(signal.SIGALRM, _h); signal.alarm(240)
                    try:
                        cert = certify(ex, Ot if isinstance(Ot, str) else [float(Ot[0]), float(Ot[1])], want=tc, verbose=False)
                    finally:
                        signal.alarm(0)
                    ok = (cert['count'] == tc)
                    log(f"  exact recheck of true={tc} < proxy={b}: certified count {cert['count']}")
            except Exception as e:
                log(f"  exact recheck of true={tc} < proxy={b} failed: {e!r}")
            if not ok:
                tc, Ot, d = b, U.label(int(allbest[0][2])) if allbest[0][2] >= 0 else None, None
                log(f"  -> using the proxy value {b} (uncertified float centre discarded)")
        prev_best = [int(i) for i in S]
        rec = dict(universe=U.name, n=n, N=U.N, M=U.M, formula=formula(n), proxy=int(b), true=int(tc),
                   nblocks=int(nb), deg=(int(d) if d is not None else None), structure=describe(U, S),
                   centre=(None if Ot is None else (Ot if isinstance(Ot, str) else [float(Ot[0]), float(Ot[1])])),
                   proxy_centre=U.label(int(allbest[0][2])) if allbest[0][2] >= 0 else None,
                   S=[U.label(int(i)) for i in S], S_idx=[int(i) for i in S], exact=exact_coords(U, S),
                   seconds=round(time.time() - t, 1), iters=iters, restarts=restarts)
        flag = "BELOW f!" if tc < formula(n) else ("= f" if tc == formula(n) else f"+{tc - formula(n)}")
        log(f"  n={n:2d} f={formula(n):3d} proxy={b:3d} true={tc:3d} [{flag}] structure={rec['structure']} "
            f"deg={d} centre={rec['centre']} t={rec['seconds']}s")
        res[n] = rec
        with open(f"{HERE}/runs/sa_results.jsonl", "a") as fh:
            fh.write(json.dumps(rec) + "\n")
        if tc <= formula(n):
            with open(f"{HERE}/runs/sa_records.jsonl", "a") as fh:
                fh.write(json.dumps(rec) + "\n")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('names', nargs='+')
    ap.add_argument('--ns', default='9-16')
    ap.add_argument('--iters', type=int, default=300_000)
    ap.add_argument('--restarts', type=int, default=6)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--closure', action='store_true', help='also run the coincidence-closed universe')
    ap.add_argument('--tlimit', type=float, default=None, help='seconds per n')
    ap.add_argument('--skip-done', action='store_true', help='skip (universe, n) pairs already in runs/sa_results.jsonl')
    ap.add_argument('--budget', type=float, default=4e8, help='block visits per restart (caps iterations for dense universes)')
    a = ap.parse_args()
    lo, hi = map(int, a.ns.split('-'))
    ns = list(range(lo, hi + 1))
    C = catalogue()
    names = list(C.keys()) if a.names == ['all'] else a.names
    for name in names:
        t = time.time()
        try:
            U = C[name]()
        except Exception as e:
            print(f"[{name}] build failed: {e!r}", flush=True)
            continue
        print(f"[{name}] {U.summary()}", flush=True)
        run_universe(U, ns, a.iters, a.restarts, a.seed, lambda s: print(s, flush=True), a.tlimit, skip_done=a.skip_done, budget=a.budget)
        if a.closure:
            try:
                U2 = closure_universe(U, name + '+X')
                print(f"[{name}+X] {U2.summary()}", flush=True)
                run_universe(U2, ns, a.iters, a.restarts, a.seed, lambda s: print(s, flush=True), a.tlimit, skip_done=a.skip_done, budget=a.budget)
            except Exception as e:
                print(f"[{name}+X] failed: {e!r}", flush=True)
        print(f"[{name}] done in {time.time()-t:.0f}s", flush=True)
