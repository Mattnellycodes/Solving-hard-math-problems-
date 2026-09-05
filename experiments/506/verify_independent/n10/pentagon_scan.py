"""Main agent's deterministic scan (2026-09-05) of the decisive geometric step for c(10) = 33.
Structure: the theorem-only survivor of the n10-continue enumeration (runs/post_7.json entry 4): two disjoint
5-blocks S, R + 20 four-blocks + 10 three-point lines, count 32. Its realisations (agent's torsion analysis):
S on the unit circle at angles 2*pi*lab/10, R at radius rho with the same angle labelling (8 labellings).
Claim to test: for no rho (rho != 0, 1) does a point O lie on the circumcircles of all 10 line-triples
(O would be the inversion centre making them straight lines). Method: for each labelling and rho on a fine
grid, intersect the first two circumcircles and measure the max distance of each intersection point to the
other eight circles; also test O = infinity (all triples collinear). Report minima away from rho = 1.
"""
import json, cmath, math, itertools
post = json.load(open('/home/user/Solving-hard-math-problems-/experiments/506/n10-continue/runs/post_7.json')); rec = post[4]
F = [frozenset(B) for B in rec['blocks']]; L = [sorted(S) for S in rec['lines']]
fives = sorted([B for B in F if len(B) == 5], key=sorted); S, R = sorted(fives[0]), sorted(fives[1])
idx = {p: i for i, p in enumerate(S + R)}
labellings = [(0,2,6,4,8,1,7,3,9,5),(0,2,6,4,8,6,2,8,4,0),(0,4,2,8,6,2,4,6,8,0),(0,4,2,8,6,7,9,1,3,5),(0,6,8,2,4,3,1,9,7,5),(0,6,8,2,4,8,6,4,2,0),(0,8,4,6,2,4,8,2,6,0),(0,8,4,6,2,9,3,7,1,5)]
print("blocks:", len(F), "lines:", L, "count:", 120 - sum(math.comb(len(B),3)-1 for B in F) - len(L))
def pts(lab, rho):
    z = {}
    for p in S: z[p] = cmath.exp(2j*math.pi*lab[idx[p]]/10)
    for p in R: z[p] = rho*cmath.exp(2j*math.pi*lab[idx[p]]/10)
    return z
def circ(a,b,c):
    ax,ay=a.real,a.imag; bx,by=b.real,b.imag; cx,cy=c.real,c.imag
    d=2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    if abs(d)<1e-13: return None
    a2=ax*ax+ay*ay; b2=bx*bx+by*by; c2=cx*cx+cy*cy
    ux=(a2*(by-cy)+b2*(cy-ay)+c2*(ay-by))/d; uy=(a2*(cx-bx)+b2*(ax-cx)+c2*(bx-ax))/d
    return (ux,uy,math.hypot(ax-ux,ay-uy))
def inter(c0,c1):
    (x0,y0,r0),(x1,y1,r1)=c0,c1; D=math.hypot(x1-x0,y1-y0)
    if D<1e-13 or D>r0+r1 or D<abs(r0-r1): return []
    a=(r0*r0-r1*r1+D*D)/(2*D); h=math.sqrt(max(r0*r0-a*a,0)); mx=x0+a*(x1-x0)/D; my=y0+a*(y1-y0)/D
    return [(mx+h*(y1-y0)/D, my-h*(x1-x0)/D),(mx-h*(y1-y0)/D, my+h*(x1-x0)/D)]
def rel_err(O, c):  # distance to circle relative to its radius (scale-free)
    return abs(math.hypot(O[0]-c[0],O[1]-c[1])-c[2])/c[2]
import numpy as np
for k, lab in enumerate(labellings):
    z = pts(lab, 0.5)
    assert all(circ(*[z[p] for p in sorted(B)[:3]]) and abs(math.hypot(z[sorted(B)[3]].real-circ(*[z[p] for p in sorted(B)[:3]])[0], z[sorted(B)[3]].imag-circ(*[z[p] for p in sorted(B)[:3]])[1]) - circ(*[z[p] for p in sorted(B)[:3]])[2]) < 1e-9 for B in F if len(B)==4), "4-blocks not concyclic"
    best = {}
    for lo, hi, name in ((0.01, 0.97, 'rho<1'), (1.03, 60.0, 'rho>1'), (0.97, 1.03, 'near 1')):
        bm = (1e9, None)
        grid = np.linspace(lo, hi, 20000) if name != 'rho>1' else np.exp(np.linspace(math.log(lo), math.log(hi), 20000))
        for rho in grid:
            z = pts(lab, rho); cs = [circ(*[z[p] for p in T]) for T in L]
            if any(c is None for c in cs): continue
            # try all pairs among the first three circles as the defining pair (robustness)
            for i, j in ((0,1),(0,2),(1,2)):
                for O in inter(cs[i], cs[j]):
                    e = max(rel_err(O, c) for t, c in enumerate(cs) if t not in (i, j))
                    if e < bm[0]: bm = (e, float(rho))
        best[name] = bm
    # O = infinity: all 10 triples collinear? measure the max normalised collinearity defect over rho
    inf_best = (1e9, None)
    for rho in np.linspace(0.01, 0.99, 5000):
        z = pts(lab, rho); e = 0
        for T in L:
            a, b, c = (z[p] for p in T); d = (b-a)*(c-a).conjugate(); e = max(e, abs(d.imag)/(abs(d)+1e-300))
        if e < inf_best[0]: inf_best = (e, float(rho))
    print(f"labelling {k}: min max-relative-distance to the other 8 circles: rho<1: {best['rho<1'][0]:.4f} (at rho={best['rho<1'][1]:.4f}); rho>1: {best['rho>1'][0]:.4f} (at {best['rho>1'][1]:.3f}); near rho=1: {best['near 1'][0]:.2e} (at {best['near 1'][1]:.4f}); O=inf collinearity defect min {inf_best[0]:.3f}")
