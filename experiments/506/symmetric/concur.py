"""Concurrency events: parameter values where three blocks pass through a common point
(these raise deg(O) and are invisible to quadruple-coincidence conditions).

For a 1-parameter slice we label blocks by their lexicographically least triple, trace the two
intersection points of every (rep block, block) pair along a grid, and detect sign changes of
the signed distance from these points to every third block; each sign change is refined by
bisection on t and the configuration is evaluated there."""
import sys, os, math, time, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import symcore as S


def block_labels(pts, has_inf, rep_trips):
    """blocks at pts; returns (labels (nb,3) triples, geoms (nb,4) [type, a,b,c], is_rep (nb,) bool).
    type 0 = circle (cx,cy,r), 1 = line (A,B,C)."""
    Bl = S.blocks_of(pts, has_inf)
    labs = []; geoms = []; isrep = []
    for key, mem in Bl.items():
        fin = sorted(m for m in mem if m >= 0)
        if len(fin) < 3: continue        # 2-point lines through infinity: not blocks of >=3 finite pts... keep? (deg counts them)
        lab = tuple(fin[:3])
        g = S.precise_geom(pts, fin)
        labs.append(lab); geoms.append([0 if g[0] == 'C' else 1, g[1], g[2], g[3]])
        isrep.append(lab in rep_trips)
    order = np.argsort([l[0] * 10**6 + l[1] * 10**3 + l[2] for l in labs])
    labs = [labs[o] for o in order]
    return labs, np.array([geoms[o] for o in order]), np.array([isrep[o] for o in order])


def pair_points(G, I, J):
    """intersection points of blocks I (indices) with blocks J, two per pair: returns (k,2,2) array
    and a validity mask (k,2).  Circle-circle, circle-line, line-line handled."""
    ti = G[I, 0]; tj = G[J, 0]
    out = np.full((len(I), 2, 2), np.nan)
    # circle-circle
    m = (ti == 0) & (tj == 0)
    if m.any():
        cx, cy, r = G[I[m], 1], G[I[m], 2], G[I[m], 3]
        dx, dy, r2 = G[J[m], 1], G[J[m], 2], G[J[m], 3]
        D = np.hypot(dx - cx, dy - cy)
        ok = (D > 1e-12) & (D <= r + r2) & (D >= np.abs(r - r2))
        a = (r * r - r2 * r2 + D * D) / (2 * np.where(D > 0, D, 1))
        h = np.sqrt(np.maximum(r * r - a * a, 0))
        mx, my = cx + a * (dx - cx) / np.where(D > 0, D, 1), cy + a * (dy - cy) / np.where(D > 0, D, 1)
        ux, uy = (dy - cy) / np.where(D > 0, D, 1), (dx - cx) / np.where(D > 0, D, 1)
        p = np.stack([np.stack([mx + h * ux, my - h * uy], -1), np.stack([mx - h * ux, my + h * uy], -1)], 1)
        p[~ok] = np.nan
        out[m] = p
    # circle-line (either order)
    for (mc, ci, li) in (((ti == 0) & (tj == 1), I, J), ((ti == 1) & (tj == 0), J, I)):
        if mc.any():
            cx, cy, r = G[ci[mc], 1], G[ci[mc], 2], G[ci[mc], 3]
            A, B, C = G[li[mc], 1], G[li[mc], 2], G[li[mc], 3]
            d = A * cx + B * cy - C
            ok = np.abs(d) <= r
            h = np.sqrt(np.maximum(r * r - d * d, 0))
            px, py = cx - A * d, cy - B * d
            p = np.stack([np.stack([px - B * h, py + A * h], -1), np.stack([px + B * h, py - A * h], -1)], 1)
            p[~ok] = np.nan
            out[mc] = p
    m = (ti == 1) & (tj == 1)
    if m.any():
        A1, B1, C1 = G[I[m], 1], G[I[m], 2], G[I[m], 3]
        A2, B2, C2 = G[J[m], 1], G[J[m], 2], G[J[m], 3]
        det = A1 * B2 - A2 * B1
        ok = np.abs(det) > 1e-12
        det = np.where(ok, det, 1)
        p = np.stack([(C1 * B2 - C2 * B1) / det, (A1 * C2 - A2 * C1) / det], -1)
        p[~ok] = np.nan
        out[m, 0] = p
    return out


def signed_dist(G, X):
    """signed distance from points X (k,2) to every block: (k, nb)."""
    t = G[:, 0]
    dc = np.hypot(X[:, None, 0] - G[None, :, 1], X[:, None, 1] - G[None, :, 2]) - G[None, :, 3]
    dl = X[:, None, 0] * G[None, :, 1] + X[:, None, 1] * G[None, :, 2] - G[None, :, 3]
    return np.where(t[None, :] == 0, dc, dl)


def slice_data(fam, params, rep_trips):
    pts = fam.points(params)
    labs, G, isrep = block_labels(pts, fam.has_inf, rep_trips)
    nb = len(labs)
    I0 = np.where(isrep)[0]
    I = np.repeat(I0, nb); J = np.tile(np.arange(nb), len(I0))
    keep = I != J
    I, J = I[keep], J[keep]
    P = pair_points(G, I, J)              # (k,2,2)
    D = np.stack([signed_dist(G, P[:, 0]), signed_dist(G, P[:, 1])], 1)   # (k,2,nb)
    # exclude the two defining blocks, and intersection points that are points of P
    D[np.arange(len(I)), :, I] = np.nan; D[np.arange(len(I)), :, J] = np.nan
    for w in (0, 1):
        dP = np.nanmin(np.hypot(P[:, w, None, 0] - pts[None, :, 0], P[:, w, None, 1] - pts[None, :, 1]), axis=1)
        D[dP < 1e-6, w, :] = np.nan
    return labs, I, J, D


def geom_from_triple(pts, lab):
    g = S.precise_geom(pts, list(lab))
    return np.array([0 if g[0] == 'C' else 1, g[1], g[2], g[3]])


def dist_fn(fam, params, l1, l2, w, l3):
    pts = fam.points(params)
    G = np.stack([geom_from_triple(pts, l1), geom_from_triple(pts, l2), geom_from_triple(pts, l3)])
    P = pair_points(G, np.array([0]), np.array([1]))[0, w]
    if np.isnan(P).any(): return np.nan
    return signed_dist(G, P[None, :])[0, 2]


def concurrency_roots(fam, base, i, G=240, log=None, max_cands=4000):
    base = np.asarray(base, float)
    grid = fam.grid(fam.kinds[i], G)
    # representative triples (finite points only) under the group
    rep_trips = set(tuple(t) for t in fam.trips if max(t) < fam.n)
    prev = None; roots = []
    ncand = 0
    for gi, t in enumerate(grid):
        p = base.copy(); p[i] = t
        try:
            cur = slice_data(fam, p, rep_trips)
        except Exception:
            prev = None; continue
        if prev is not None and prev[0] == cur[0] and prev[3].shape == cur[3].shape:
            D0, D1 = prev[3], cur[3]
            s = np.sign(D0) * np.sign(D1)
            small = (np.abs(D0) + np.abs(D1)) < 0.05
            notzero = np.maximum(np.abs(D0), np.abs(D1)) > 1e-9
            cand = np.argwhere((s < 0) & small & notzero & np.isfinite(D0) & np.isfinite(D1))
            # dedupe by interpolated root (symmetric copies of one event give the same t*)
            groups = {}
            for (k, w, b3) in cand:
                d0, d1 = D0[k, w, b3], D1[k, w, b3]
                tt = grid[gi - 1] + (grid[gi] - grid[gi - 1]) * d0 / (d0 - d1)
                key = round(tt / (grid[gi] - grid[gi - 1]) * 200)
                groups.setdefault(key, (k, w, b3))
            for key, (k, w, b3) in list(groups.items())[:max_cands]:
                l1, l2, l3 = cur[0][cur[1][k]], cur[0][cur[2][k]], cur[0][b3]
                lo, hi = grid[gi - 1], grid[gi]
                flo = D0[k, w, b3]
                ok = True
                for _ in range(50):
                    mid = 0.5 * (lo + hi); p2 = base.copy(); p2[i] = mid
                    fm = dist_fn(fam, p2, l1, l2, w, l3)
                    if not np.isfinite(fm): ok = False; break
                    if np.sign(fm) == np.sign(flo) or fm == 0: lo, flo = mid, fm
                    else: hi = mid
                    if hi - lo < 1e-12: break
                if ok and abs(flo) < 1e-6:
                    roots.append(0.5 * (lo + hi)); ncand += 1
        prev = cur
    roots = sorted(set(round(float(r), 12) for r in roots))
    if log: log(f"    concurrency p{i}: {ncand} candidates -> {len(roots)} distinct roots")
    return roots


if __name__ == '__main__':
    from search import Family
    fam = Family(2, [('G',), ('G',)])
    t0 = time.time()
    r = concurrency_roots(fam, np.array([1.0, math.pi / 4, 3.0, math.pi / 4]), 2, G=200, log=print)
    print(time.time() - t0, 's; roots near 7.873?', [x for x in r if abs(x - (4 + math.sqrt(15))) < 1e-3], 'near 0.127?', [x for x in r if abs(x - (4 - math.sqrt(15))) < 1e-3], len(r))
    for x in r:
        p = np.array([1.0, math.pi / 4, x, math.pi / 4])
        ev = S.evaluate_fast(fam.points(p), fam.has_inf)
        if ev and ev['best'] <= 18: print(round(x, 6), ev['best'], ev['deg'], ev['nblocks'])
