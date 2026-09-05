#!/usr/bin/env python3
"""c2_degree_scan.py -- numerical scan of the concentric-pentagon family realising structure C-2:
points  e^{2 pi i j/5} (j = 0..4)  and  r e^{i(2 pi k/5 + psi)} (k = 0..4),  psi in {0, pi/5},  0 < r < 1.
For each r: compute all blocks (circles/lines through >= 3 points), check that the block structure is exactly C-2
(42 blocks: 2 five-point, 20 four-point, 20 three-point), and compute the maximum number of blocks through a
common point O not in P (candidates: pairwise intersections of blocks).  The planar circle count after
inversion in O is 42 - deg(O).  Isolated special values of r are located by tracking, for every (pair of blocks,
third block), the distance from the intersection point to the third block and refining its local minima."""
import sys, math, cmath, itertools
import numpy as np

def points(r, psi):
    return [cmath.exp(2j * math.pi * j / 5) for j in range(5)] + [r * cmath.exp(1j * (2 * math.pi * k / 5 + psi)) for k in range(5)]

def circle3(a, b, c):
    d = 2 * (a.real * (b.imag - c.imag) + b.real * (c.imag - a.imag) + c.real * (a.imag - b.imag))
    if abs(d) < 1e-12:
        return None
    a2, b2, c2 = abs(a) ** 2, abs(b) ** 2, abs(c) ** 2
    ux = (a2 * (b.imag - c.imag) + b2 * (c.imag - a.imag) + c2 * (a.imag - b.imag)) / d
    uy = (a2 * (c.real - b.real) + b2 * (a.real - c.real) + c2 * (b.real - a.real)) / d
    ctr = complex(ux, uy)
    return ctr, abs(a - ctr)

def blocks(pts, tol=1e-9):
    """list of (kind, params, set of point indices); kind 'C' (centre, radius) or 'L' (unit normal, offset)."""
    out = []
    n = len(pts)
    for i, j, k in itertools.combinations(range(n), 3):
        if any(i in S and j in S and k in S for _, _, S in out):
            continue
        cc = circle3(pts[i], pts[j], pts[k])
        if cc is None:
            dvec = (pts[j] - pts[i]); nrm = complex(-dvec.imag, dvec.real) / abs(dvec)
            off = (nrm.conjugate() * pts[i]).real
            S = {q for q in range(n) if abs((nrm.conjugate() * pts[q]).real - off) < tol}
            out.append(('L', (nrm, off), S))
        else:
            ctr, rad = cc
            S = {q for q in range(n) if abs(abs(pts[q] - ctr) - rad) < tol}
            out.append(('C', (ctr, rad), S))
    return out

def dist_to_block(O, blk):
    kind, par, _ = blk
    if kind == 'C':
        return abs(abs(O - par[0]) - par[1])
    nrm, off = par
    return abs((nrm.conjugate() * O).real - off)

def intersections(b1, b2):
    k1, p1, _ = b1; k2, p2, _ = b2
    if k1 == 'C' and k2 == 'C':
        c1, r1 = p1; c2, r2 = p2
        d = abs(c2 - c1)
        if d < 1e-12 or d > r1 + r2 + 1e-9 or d < abs(r1 - r2) - 1e-9:
            return []
        aa = (r1 * r1 - r2 * r2 + d * d) / (2 * d); h = math.sqrt(max(r1 * r1 - aa * aa, 0.0))
        p = c1 + aa * (c2 - c1) / d
        return [p + h * 1j * (c2 - c1) / d, p - h * 1j * (c2 - c1) / d]
    if k1 == 'L' and k2 == 'L':
        n1, o1 = p1; n2, o2 = p2
        det = n1.real * n2.imag - n1.imag * n2.real
        if abs(det) < 1e-12: return []
        x = (o1 * n2.imag - o2 * n1.imag) / det; y = (n1.real * o2 - n2.real * o1) / det
        return [complex(x, y)]
    if k1 == 'L': b1, b2 = b2, b1; k1, p1, _ = b1; k2, p2, _ = b2
    c, rad = p1; nrm, off = p2
    # foot of perpendicular from c to the line
    dist = (nrm.conjugate() * c).real - off
    foot = c - dist * nrm
    h2 = rad * rad - dist * dist
    if h2 < -1e-9: return []
    h = math.sqrt(max(h2, 0.0)); t = complex(-nrm.imag, nrm.real)
    return [foot + h * t, foot - h * t]

def analyse(r, psi, tol=1e-7):
    pts = points(r, psi)
    bl = blocks(pts)
    sizes = sorted(len(S) for _, _, S in bl)
    exact = (len(bl) == 42 and sizes.count(5) == 2 and sizes.count(4) == 20 and sizes.count(3) == 20)
    best = (0, None)
    for b1, b2 in itertools.combinations(bl, 2):
        for O in intersections(b1, b2):
            if min(abs(O - q) for q in pts) < 1e-6: continue
            deg = sum(1 for b in bl if dist_to_block(O, b) < tol)
            if deg > best[0]: best = (deg, O)
    nlines = sum(1 for k, _, _ in bl if k == 'L')
    return exact, len(bl), sizes, best, nlines

if __name__ == "__main__":
    for psi_name, psi in (("0", 0.0), ("pi/5", math.pi / 5)):
        print(f"=== psi = {psi_name}")
        grid = np.linspace(0.02, 0.98, 481)
        summary = {}
        for r in grid:
            exact, nb, sizes, (deg, O), nl = analyse(float(r), psi)
            key = (exact, nb, deg)
            summary.setdefault(key, []).append(round(float(r), 4))
        for key, rs in sorted(summary.items()):
            print(f"   exact-C2={key[0]} #blocks={key[1]} max deg(O)={key[2]}: {len(rs)} grid values, e.g. {rs[:5]}")
        # locate special r: minimise, over r, the distance from intersection points of block pairs to third blocks.
        # We do this by a fine scan of the 'max deg' function is insufficient for isolated r, so track the
        # function g(r) = min over (pair, third, +-) of the distance, per triple, on a fine grid + refinement.
        fine = np.linspace(0.02, 0.98, 4801)
        special = {}
        prev = None
        print("   scanning for isolated special r (triple concurrency) ...", flush=True)
        # per r compute the list of (pair index, third index, sign) distances; keep per-triple minima across the grid
        trip_min = {}
        for r in fine:
            pts = points(float(r), psi); bl = blocks(pts)
            if len(bl) != 42: 
                special[round(float(r), 6)] = ('extra incidences', len(bl)); continue
            # blocks are generated in a canonical order (by first triple) -> stable indexing across r
            for i, j in itertools.combinations(range(42), 2):
                for s, O in enumerate(intersections(bl[i], bl[j])):
                    if min(abs(O - q) for q in pts) < 1e-6: continue
                    for k in range(42):
                        if k in (i, j): continue
                        d = dist_to_block(O, bl[k])
                        key = (i, j, s, k)
                        if key not in trip_min or d < trip_min[key][0]:
                            trip_min[key] = (d, float(r))
        cands = sorted(set(round(v[1], 6) for v in trip_min.values() if v[0] < 1e-3))
        print(f"   candidate special r values (triple nearly concurrent on the fine grid): {len(cands)}")
        # refine each candidate: golden-section search of max-degree count is discrete; instead refine the specific
        # triple distance by bisection-like minimisation, then compute the degree at the refined r.
        import scipy.optimize as so
        refined = {}
        for key, (d, r0) in trip_min.items():
            if d >= 1e-3: continue
            i, j, s, k = key
            def f(r):
                pts = points(float(r), psi); bl = blocks(pts)
                if len(bl) != 42: return 1.0
                I = intersections(bl[i], bl[j])
                if len(I) <= s: return 1.0
                return dist_to_block(I[s], bl[k])
            res = so.minimize_scalar(f, bounds=(max(0.02, r0 - 0.001), min(0.98, r0 + 0.001)), method='bounded',
                                     options={'xatol': 1e-12})
            if res.fun < 1e-8:
                rr = round(float(res.x), 9)
                refined.setdefault(rr, []).append(key)
        print(f"   refined special r values: {sorted(refined)}")
        for rr in sorted(refined):
            exact, nb, sizes, (deg, O), nl = analyse(rr, psi, tol=1e-6)
            print(f"      r = {rr:.9f}: exact-C2={exact} #blocks={nb} sizes={sorted(set(sizes))} max deg(O)={deg} at O={O}  (#lines={nl})")
