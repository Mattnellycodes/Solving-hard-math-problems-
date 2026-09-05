"""High-precision (mpmath, 40 digits) block structure of two concentric regular pentagons (same orientation)."""
import mpmath as mp
from itertools import combinations
from math import comb
mp.mp.dps = 40
def blocks(pts, tol=mp.mpf('1e-25')):
    n = len(pts); out = {}
    for t in combinations(range(n), 3):
        (ax, ay), (bx, by), (cx, cy) = (pts[i] for i in t)
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        if abs(d) < tol:
            key = ('L',) + tuple(sorted(t))  # generic: no 4 collinear here; verify below
            blk = set(t)
            for m in range(n):
                if m in blk: continue
                mx, my = pts[m]
                if abs((bx-ax)*(my-ay)-(by-ay)*(mx-ax)) < tol: blk.add(m)
            out[('L', frozenset(blk))] = 1
        else:
            a2, b2, c2 = ax*ax+ay*ay, bx*bx+by*by, cx*cx+cy*cy
            ux = (a2*(by-cy)+b2*(cy-ay)+c2*(ay-by))/d; uy = (a2*(cx-bx)+b2*(ax-cx)+c2*(bx-ax))/d
            r2 = (ax-ux)**2+(ay-uy)**2
            blk = set(t)
            for m in range(n):
                if m in blk: continue
                mx, my = pts[m]
                if abs((mx-ux)**2+(my-uy)**2 - r2) < tol * (1 + r2): blk.add(m)
            out[('C', frozenset(blk))] = 1
    return [k for k in out]
for r in (mp.mpf('1.7320508'), mp.mpf('2.3'), mp.mpf('0.61803398874989484820458683436563811772')):
    pts = [(mp.cos(2*mp.pi*k/5), mp.sin(2*mp.pi*k/5)) for k in range(5)] + [(r*mp.cos(2*mp.pi*k/5), r*mp.sin(2*mp.pi*k/5)) for k in range(5)]
    B = blocks(pts)
    circ = [b for b in B if b[0] == 'C']; lin = [b for b in B if b[0] == 'L']
    D = sum(comb(len(b[1]), 3) - 1 for b in B if len(b[1]) >= 4)
    from collections import Counter
    print(f"r={mp.nstr(r,8)}: circles={len(circ)} lines={len(lin)} D={D} sizes={Counter(len(b[1]) for b in B)}  check C(10,3)-D-l={comb(10,3)-D-len(lin)}")
