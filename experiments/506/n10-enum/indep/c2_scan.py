"""Numeric scan (evidence): for each realisation type of C-2 (from c2_rigidity.json) and rho on a grid,
compute the 42 blocks (22 rich blocks + 20 uncovered triples) as circles/lines, all pairwise
intersection points, and the maximum number of blocks through a common point O not in P."""
import json, math, itertools, cmath
import numpy as np
import mo
C2 = json.load(open('c2_rigidity.json'))
blocks = [tuple(b) for b in C2['blocks']]
F = [mo.mask(b) for b in blocks]
allblocks = blocks + [tuple(mo.bits(t)) for t in mo.uncovered_triples(10, F)]
assert len(allblocks) == 42
A, B = C2['A'], C2['B']

def circle_through(p, q, r):
    """generalized circle a(x^2+y^2) + b x + c y + d = 0 through three complex points (normalised)."""
    M = np.array([[abs(z) ** 2, z.real, z.imag, 1.0] for z in (p, q, r)])
    # null vector of the 3x4 matrix
    u, s, vt = np.linalg.svd(M)
    v = vt[-1]
    return v / np.linalg.norm(v)

def intersections(c1, c2):
    a1, b1, c1_, d1 = c1; a2, b2, c2_, d2 = c2
    # eliminate the quadratic term: line L = a2*c1 - a1*c2 (radical axis), then intersect with circle 1 (or 2)
    L = a2 * c1 - a1 * c2   # coefficients (0, bx, cy, d)
    _, lb, lc, ld = L
    if abs(lb) < 1e-12 and abs(lc) < 1e-12:
        return []
    # parametrise the line
    if abs(lc) >= abs(lb):
        # y = -(lb x + ld)/lc
        pts = []
        # substitute in circle with larger |a|
        a, b, c, d = c1 if abs(a1) >= abs(a2) else c2
        # x^2 + y^2 terms: a(x^2 + y^2) + b x + c y + d = 0 with y = m x + k
        m = -lb / lc; k = -ld / lc
        qa = a * (1 + m * m); qb = 2 * a * m * k + b + c * m; qc = a * k * k + c * k + d
        if abs(qa) < 1e-14:
            if abs(qb) < 1e-14: return []
            x = -qc / qb; return [complex(x, m * x + k)]
        disc = qb * qb - 4 * qa * qc
        if disc < -1e-9: return []
        disc = max(disc, 0.0)
        for sgn in (1, -1):
            x = (-qb + sgn * math.sqrt(disc)) / (2 * qa)
            pts.append(complex(x, m * x + k))
        return pts
    else:
        m = -lc / lb; k = -ld / lb   # x = m y + k
        a, b, c, d = c1 if abs(a1) >= abs(a2) else c2
        qa = a * (1 + m * m); qb = 2 * a * m * k + c + b * m; qc = a * k * k + b * k + d
        if abs(qa) < 1e-14:
            if abs(qb) < 1e-14: return []
            y = -qc / qb; return [complex(m * y + k, y)]
        disc = qb * qb - 4 * qa * qc
        if disc < -1e-9: return []
        disc = max(disc, 0.0)
        return [complex(m * y + k, y) for y in [(-qb + s * math.sqrt(disc)) / (2 * qa) for s in (1, -1)]]

def on_circle(c, z, tol=1e-7):
    a, b, cc, d = c
    return abs(a * abs(z) ** 2 + b * z.real + cc * z.imag + d) < tol * (1 + abs(z) ** 2)

results = {}
for ti, T in enumerate(C2['types']):
    ang = {p: math.radians(a) for p, a in zip(A, T['anglesA'])}
    ang.update({p: math.radians(b) for p, b in zip(B, T['anglesB'])})
    hist = {}
    worst = (0, None, None)
    for rho in np.concatenate([np.linspace(1.002, 1.5, 400), np.linspace(1.5, 6, 400), np.linspace(6, 60, 200)]):
        P = {p: cmath.exp(1j * ang[p]) for p in A}
        P.update({p: rho * cmath.exp(1j * ang[p]) for p in B})
        circ = [circle_through(P[b[0]], P[b[1]], P[b[2]]) for b in allblocks]
        # sanity: all points of each block on its circle
        assert all(on_circle(circ[i], P[q]) for i, b in enumerate(allblocks) for q in b)
        best = 0; bestpt = None
        cand_pts = []
        for i, j in itertools.combinations(range(42), 2):
            for z in intersections(circ[i], circ[j]):
                if any(abs(z - P[q]) < 1e-6 for q in range(10)):
                    continue
                cand_pts.append(z)
        for z in cand_pts:
            k = sum(1 for c in circ if on_circle(c, z, 1e-6))
            if k > best:
                best = k; bestpt = z
        hist[best] = hist.get(best, 0) + 1
        if best > worst[0]:
            worst = (best, float(rho), bestpt)
    print(f"type {ti} (B offset {T['B_offset_deg']} deg): histogram of max #blocks through a non-configuration point over rho grid: {dict(sorted(hist.items()))}; max = {worst[0]} at rho = {worst[1]}")
    results[ti] = {'hist': {int(k): v for k, v in hist.items()}, 'max': worst[0], 'rho': worst[1]}
json.dump(results, open('c2_scan.json', 'w'), indent=1)
