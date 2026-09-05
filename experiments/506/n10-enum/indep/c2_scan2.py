"""Robust numeric scan (Moebius-normalised tests): max number of the 42 blocks of C-2 through a point
O not in P, for both geometric types and rho on a grid.  Evidence only."""
import json, math, itertools, cmath
import numpy as np
import mo
rig = json.load(open('c2_rigidity.json'))
blocks = [tuple(b) for b in rig['blocks']]
F = [mo.mask(b) for b in blocks]
allblocks = blocks + [tuple(mo.bits(t)) for t in mo.uncovered_triples(10, F)]
A, B = rig['A'], rig['B']
types = {}
for T in rig['types']:
    types.setdefault(T['B_offset_deg'], T)

def mob(p, q, r):
    k = (q - r) / (q - p)
    fwd = lambda z: k * (z - p) / (z - r)
    inv = lambda w: (w * r - k * p) / (w - k)
    return fwd, inv

def on_block(P, bl, z, tol=1e-9):
    p, q, r = (P[i] for i in bl[:3])
    if min(abs(z - p), abs(z - q), abs(z - r)) < 1e-12:
        return True
    cr = ((z - p) * (q - r)) / ((z - r) * (q - p))
    return abs(cr.imag) < tol * (1 + abs(cr))

def inter(P, bi, bj):
    shared = [i for i in bi if i in bj]
    order = shared + [i for i in bi if i not in bj]
    p, q, r = (P[i] for i in order[:3])
    fwd, inv = mob(p, q, r)
    Pp, Qp, Rp = (fwd(P[i]) for i in bj[:3])
    u = (Qp - Rp) / (Qp - Pp)
    ca = u.imag; cb = -(u * (Pp + Rp.conjugate())).imag; cc = (u * Pp * Rp.conjugate()).imag
    roots = []
    if abs(ca) < 1e-14:
        if abs(cb) > 1e-14:
            roots = [-cc / cb]
    else:
        disc = cb * cb - 4 * ca * cc
        if disc >= 0:
            roots = [(-cb + sg * math.sqrt(disc)) / (2 * ca) for sg in (1, -1)]
    out = []
    for w in roots:
        try:
            out.append(inv(complex(w, 0)))
        except ZeroDivisionError:
            pass
    return out

res = {}
for off, T in types.items():
    ang = {p: math.radians(a) for p, a in zip(A, T['anglesA'])}
    ang.update({p: math.radians(b) for p, b in zip(B, T['anglesB'])})
    hist = {}; worst = (0, None); special = []
    grid = np.concatenate([np.linspace(0.02, 0.98, 300), np.linspace(1.02, 3, 400), np.linspace(3, 12, 300), np.linspace(12, 60, 100)])
    for r_ in grid:
        P = {p: cmath.exp(1j * ang[p]) for p in A}
        P.update({p: r_ * cmath.exp(1j * ang[p]) for p in B})
        pts = []
        for i, j in itertools.combinations(range(42), 2):
            for z in inter(P, allblocks[i], allblocks[j]):
                if all(abs(z - P[q]) > 1e-7 for q in range(10)):
                    pts.append(z)
        best = 0
        for z in pts:
            k = sum(1 for bl in allblocks if on_block(P, bl, z))
            best = max(best, k)
        hist[best] = hist.get(best, 0) + 1
        if best >= 5:
            special.append((round(float(r_), 4), best))
        if best > worst[0]:
            worst = (best, float(r_))
    print(f"type offset {off}: histogram of max #blocks through a non-configuration point: {dict(sorted(hist.items()))}; max {worst[0]} at rho={worst[1]}; grid values with >=5: {special[:20]}", flush=True)
    res[off] = {'hist': {int(k): v for k, v in hist.items()}, 'max': worst[0], 'rho': worst[1], 'special': special}
json.dump(res, open('c2_scan2.json', 'w'), indent=1)
