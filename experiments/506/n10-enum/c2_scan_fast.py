#!/usr/bin/env python3
"""c2_scan_fast.py -- vectorised numerical scan of the concentric-pentagon family realising structure C-2:
points e^{2 pi i j/5} (j = 0..4) and r e^{i(2 pi k/5 + psi)} (k = 0..4), psi in {0, pi/5}, 0 < r < 1.
For each r on a fine grid we compute the 42 blocks (as circles; we verify that the structure is exactly C-2),
all pairwise intersection points, and the maximum number of blocks through one intersection point.
To catch ISOLATED special values of r (where an extra block passes through an intersection point), we track for
every (ordered pair of blocks, intersection branch, third block) the distance function along the grid and refine
every local minimum below a threshold by bounded scalar minimisation; at each refined r we recompute the maximum
degree with a tight tolerance.  Output: the maximum deg(O) over the whole family (the planar circle count of a
realisation is 42 - deg(O)) and the list of special r values with their degrees."""
import math, itertools, sys
import numpy as np
from scipy.optimize import minimize_scalar

ZETA = np.exp(2j * np.pi * np.arange(5) / 5)


def points(r, psi):
    return np.concatenate([ZETA, r * ZETA * np.exp(1j * psi)])


TRIPLES = np.array(list(itertools.combinations(range(10), 3)))


def circles_all(pts):
    """circumcircles of all 120 triples: centres (complex) and radii; NaN for collinear."""
    a, b, c = pts[TRIPLES[:, 0]], pts[TRIPLES[:, 1]], pts[TRIPLES[:, 2]]
    d = 2 * (a.real * (b.imag - c.imag) + b.real * (c.imag - a.imag) + c.real * (a.imag - b.imag))
    with np.errstate(divide='ignore', invalid='ignore'):
        a2, b2, c2 = abs(a) ** 2, abs(b) ** 2, abs(c) ** 2
        ux = (a2 * (b.imag - c.imag) + b2 * (c.imag - a.imag) + c2 * (a.imag - b.imag)) / d
        uy = (a2 * (c.real - b.real) + b2 * (a.real - c.real) + c2 * (b.real - a.real)) / d
    ctr = ux + 1j * uy
    rad = abs(a - ctr)
    return ctr, rad, np.abs(d) < 1e-11


def blocks(pts, tol=1e-8):
    """distinct blocks as (centre, radius, member set); collinear triples reported separately."""
    ctr, rad, lin = circles_all(pts)
    out = []
    used = np.zeros(len(TRIPLES), dtype=bool)
    lines = []
    for t in range(len(TRIPLES)):
        if used[t]:
            continue
        if lin[t]:
            lines.append(tuple(TRIPLES[t])); used[t] = True; continue
        members = {q for q in range(10) if abs(abs(pts[q] - ctr[t]) - rad[t]) < tol}
        for t2 in range(t, len(TRIPLES)):
            if set(TRIPLES[t2]) <= members:
                used[t2] = True
        out.append((ctr[t], rad[t], frozenset(members)))
    return out, lines


def intersections(C1, R1, C2, R2):
    d = abs(C2 - C1)
    if d < 1e-12 or d > R1 + R2 + 1e-9 or d < abs(R1 - R2) - 1e-9:
        return []
    aa = (R1 * R1 - R2 * R2 + d * d) / (2 * d)
    h = math.sqrt(max(R1 * R1 - aa * aa, 0.0))
    p = C1 + aa * (C2 - C1) / d
    return [p + h * 1j * (C2 - C1) / d, p - h * 1j * (C2 - C1) / d]


def analyse(r, psi, tol=1e-7):
    pts = points(r, psi)
    bl, lines = blocks(pts)
    sizes = sorted(len(S) for _, _, S in bl)
    exact = (not lines) and len(bl) == 42 and sizes.count(5) == 2 and sizes.count(4) == 20 and sizes.count(3) == 20
    ctrs = np.array([c for c, _, _ in bl]); rads = np.array([rr for _, rr, _ in bl])
    # all pairwise intersection points
    O = []
    for i, j in itertools.combinations(range(len(bl)), 2):
        O.extend(intersections(ctrs[i], rads[i], ctrs[j], rads[j]))
    if not O:
        return exact, len(bl), len(lines), 0, None
    O = np.array(O)
    # exclude points of P
    far = np.min(np.abs(O[:, None] - pts[None, :]), axis=1) > 1e-6
    O = O[far]
    dist = np.abs(np.abs(O[:, None] - ctrs[None, :]) - rads[None, :])      # (#O, #blocks)
    deg = (dist < tol).sum(axis=1)
    k = int(np.argmax(deg))
    return exact, len(bl), len(lines), int(deg[k]), O[k]


if __name__ == "__main__":
    overall = 0
    for psi_name, psi in (("0", 0.0), ("pi/5", math.pi / 5)):
        print(f"=== psi = {psi_name}", flush=True)
        grid = np.linspace(0.01, 0.99, 981)
        summary = {}
        for r in grid:
            exact, nb, nl, deg, O = analyse(float(r), psi)
            summary.setdefault((exact, nb, nl, deg), []).append(round(float(r), 3))
        for key, rs in sorted(summary.items()):
            print(f"   exact-C2={key[0]} #blocks={key[1]} #lines={key[2]} max deg(O)={key[3]}: {len(rs)} grid values, e.g. {rs[:4]}")
        # isolated special r: track distances of intersection points to third blocks along a fine grid
        fine = np.linspace(0.01, 0.99, 9801)
        prev_key = None
        best_local = {}
        for r in fine:
            pts = points(float(r), psi)
            bl, lines = blocks(pts)
            if lines or len(bl) != 42:
                best_local.setdefault(('extra', round(float(r), 6)), (0.0, float(r)))
                continue
            # canonical ordering of blocks by member set for stable indexing
            bl = sorted(bl, key=lambda x: tuple(sorted(x[2])))
            ctrs = np.array([c for c, _, _ in bl]); rads = np.array([rr for _, rr, _ in bl])
            for i, j in itertools.combinations(range(42), 2):
                I = intersections(ctrs[i], rads[i], ctrs[j], rads[j])
                for s, Op in enumerate(I):
                    if np.min(np.abs(Op - pts)) < 1e-6:
                        continue
                    d = np.abs(np.abs(Op - ctrs) - rads)
                    d[[i, j]] = 1.0
                    ks = np.flatnonzero(d < 2e-3)
                    for k in ks:
                        key = (i, j, s, int(k))
                        if key not in best_local or d[k] < best_local[key][0]:
                            best_local[key] = (float(d[k]), float(r))
        specials = {}
        for key, (d, r0) in best_local.items():
            if key[0] == 'extra':
                specials.setdefault(round(r0, 6), set()).add('extra-incidence grid point')
                continue
            i, j, s, k = key

            def f(rv):
                pts = points(float(rv), psi)
                bl, lines = blocks(pts)
                if lines or len(bl) != 42:
                    return 1.0
                bl = sorted(bl, key=lambda x: tuple(sorted(x[2])))
                I = intersections(bl[i][0], bl[i][1], bl[j][0], bl[j][1])
                if len(I) <= s:
                    return 1.0
                return abs(abs(I[s] - bl[k][0]) - bl[k][1])
            res = minimize_scalar(f, bounds=(max(0.01, r0 - 2e-4), min(0.99, r0 + 2e-4)), method='bounded',
                                  options={'xatol': 1e-13})
            if res.fun < 1e-9:
                specials.setdefault(round(float(res.x), 9), set()).add(key)
        print(f"   special r values found: {len(specials)}")
        for rr in sorted(specials):
            exact, nb, nl, deg, O = analyse(rr, psi, tol=1e-6)
            overall = max(overall, deg if exact else 0)
            print(f"      r = {rr:.9f}: exact-C2={exact} #blocks={nb} #lines={nl} max deg(O)={deg}"
                  + (f" at O = {O:.6f}" if O is not None else ""))
    print(f"MAXIMUM deg(O) over the family with exactly the C-2 structure: {overall}  -> circle count >= {42 - overall}")
