"""Two concentric regular pentagons (same orientation), 10 points: block structure at several radius
ratios, using mpmath at 50 digits with a 1e-30 tolerance (a float check of the radical-axis argument:
A_i A_j B_k B_l concyclic  <=>  A_iA_j || B_kB_l, giving 20 four-circles + two 5-circles, D = 78)."""
from mpmath import mp, mpf, cos, sin, pi, matrix, det
from itertools import combinations
from math import comb
mp.dps = 50
def pts(rho, theta=0):
    A = [(cos(2*pi*i/5), sin(2*pi*i/5)) for i in range(5)]
    B = [(rho*cos(2*pi*i/5+theta), rho*sin(2*pi*i/5+theta)) for i in range(5)]
    return A + B
def gcirc(p, q, r):  # (a,b,c,d) up to scale
    rows = [(x*x+y*y, x, y, mpf(1)) for x, y in (p, q, r)]
    def minor(cols):
        return det(matrix([[row[c] for c in cols] for row in rows]))
    return (minor((1,2,3)), -minor((0,2,3)), minor((0,1,3)), -minor((0,1,2)))
def blocks(P, tol=mpf('1e-30')):
    n = len(P); used = set(); out = []
    for i, j, k in combinations(range(n), 3):
        if (i, j, k) in used: continue
        a, b, c, d = gcirc(P[i], P[j], P[k])
        s = max(abs(a), abs(b), abs(c), abs(d))
        mem = [m for m in range(n) if abs(a*(P[m][0]**2+P[m][1]**2)+b*P[m][0]+c*P[m][1]+d)/s < tol]
        for t in combinations(mem, 3): used.add(t)
        out.append((len(mem), abs(a)/s < tol))
    assert sum(comb(k, 3) for k, _ in out) == comb(n, 3)
    return out
for rho in (mpf('0.5'), mpf('0.61803398874989484820458683436563811772'), mpf('0.73'), mpf('0.381966011250105151795413165634361882')):
    for theta in (0, pi/5):
        bl = blocks(pts(rho, theta))
        D = sum(comb(k, 3)-1 for k, _ in bl if k >= 4); l = sum(1 for _, isline in bl if isline)
        sizes = sorted((k for k, _ in bl), reverse=True)
        print(f"rho={mp.nstr(rho,8)} theta={'0' if theta==0 else 'pi/5'}: blocks={len(bl)} sizes>=4: {[k for k in sizes if k>=4]} D={D} lines={l} circles={comb(10,3)-D-l}")
