"""Numerical cross-check for the concentric-pentagon family (own code): P = {zeta^k} ∪ {rho * e^{i psi} zeta^k},
psi ∈ {0, pi/5}.  For rho on a fine grid compute all blocks (circles/lines through >= 3 points) and the maximum
number of blocks through a common point O ∉ P (O = pairwise intersections of blocks, and O = infinity); the planar
circle count after inversion in O is #blocks - deg(O).  Prints the best (smallest) counts found and the rho values
where deg(O) >= 3 occurs."""
import math, cmath, itertools, sys
import numpy as np

def points(rho, psi):
    return [cmath.exp(2j * math.pi * k / 5) for k in range(5)] + [rho * cmath.exp(1j * (2 * math.pi * k / 5 + psi)) for k in range(5)]

def blocks(pts, tol=1e-9):
    out = []; n = len(pts)
    for i, j, k in itertools.combinations(range(n), 3):
        if any({i, j, k} <= S for _, _, S in out): continue
        a, b, c = pts[i], pts[j], pts[k]
        d = 2 * (a.real * (b.imag - c.imag) + b.real * (c.imag - a.imag) + c.real * (a.imag - b.imag))
        if abs(d) < 1e-12:
            dv = b - a; nrm = complex(-dv.imag, dv.real) / abs(dv); off = (nrm.conjugate() * a).real
            S = {q for q in range(n) if abs((nrm.conjugate() * pts[q]).real - off) < tol}
            out.append(('L', (nrm, off), S))
        else:
            a2, b2, c2 = abs(a) ** 2, abs(b) ** 2, abs(c) ** 2
            ux = (a2 * (b.imag - c.imag) + b2 * (c.imag - a.imag) + c2 * (a.imag - b.imag)) / d
            uy = (a2 * (c.real - b.real) + b2 * (a.real - c.real) + c2 * (b.real - a.real)) / d
            ctr = complex(ux, uy); r = abs(a - ctr)
            S = {q for q in range(n) if abs(abs(pts[q] - ctr) - r) < tol}
            out.append(('C', (ctr, r), S))
    return out

def dist(O, blk):
    kind, par, _ = blk
    if kind == 'C': return abs(abs(O - par[0]) - par[1])
    nrm, off = par; return abs((nrm.conjugate() * O).real - off)

def inters(b1, b2):
    k1, p1, _ = b1; k2, p2, _ = b2
    if k1 == 'C' and k2 == 'C':
        c1, r1 = p1; c2, r2 = p2; d = abs(c2 - c1)
        if d < 1e-12 or d > r1 + r2 + 1e-9 or d < abs(r1 - r2) - 1e-9: return []
        aa = (r1 * r1 - r2 * r2 + d * d) / (2 * d); h = math.sqrt(max(r1 * r1 - aa * aa, 0.0))
        p = c1 + aa * (c2 - c1) / d
        return [p + h * 1j * (c2 - c1) / d, p - h * 1j * (c2 - c1) / d]
    if k1 == 'L' and k2 == 'L':
        n1, o1 = p1; n2, o2 = p2; det = n1.real * n2.imag - n1.imag * n2.real
        if abs(det) < 1e-12: return []
        return [complex((o1 * n2.imag - o2 * n1.imag) / det, (n1.real * o2 - n2.real * o1) / det)]
    if k1 == 'L': b1, b2 = b2, b1; k1, p1, _ = b1; k2, p2, _ = b2
    c, r = p1; nrm, off = p2
    dd = (nrm.conjugate() * c).real - off; foot = c - dd * nrm; h2 = r * r - dd * dd
    if h2 < -1e-9: return []
    h = math.sqrt(max(h2, 0.0)); t = complex(-nrm.imag, nrm.real)
    return [foot + h * t, foot - h * t]

def analyse(rho, psi, tol=1e-7):
    pts = points(rho, psi); bl = blocks(pts)
    if any(len(S) == 10 for _, _, S in bl): return None
    best = sum(1 for k, _, _ in bl if k == 'L'); bestO = 'inf'
    for b1, b2 in itertools.combinations(bl, 2):
        for O in inters(b1, b2):
            if any(abs(O - p) < 1e-6 for p in pts): continue
            deg = sum(1 for b in bl if dist(O, b) < tol)
            if deg > best: best, bestO = deg, O
    return len(bl), best, bestO, sorted(len(S) for _, _, S in bl)

if __name__ == '__main__':
    grid = np.concatenate([np.linspace(0.02, 0.98, 481)])
    for psi in (0.0, math.pi / 5):
        results = []
        for rho in grid:
            r = analyse(rho, psi)
            if r is None: continue
            nb, deg, O, sizes = r
            results.append((nb - deg, deg, nb, rho, O))
        results.sort(key=lambda t: (t[0], -t[1]))
        print(f"psi = {psi:.4f}: grid of {len(grid)} rho values; best (count, deg, #blocks, rho, O):")
        for t in results[:6]: print("   ", t)
        degs = sorted(set(t[1] for t in results)); print("   degrees observed:", degs, " #blocks observed:", sorted(set(t[2] for t in results)))
        hi = [t for t in results if t[1] >= 3]
        print(f"   rho values with deg(O) >= 3: {len(hi)} e.g. {[(round(t[3],4), t[1], t[2]) for t in hi[:10]]}")
