"""Skeptic c7-audit: test the Moebius identity  circles = C(n,3) - D - ell  (D over ALL blocks of
size >= 4, lines included; ell = #lines with >= 3 points) on every 7-subset of a 4x4 grid and on
random small-integer sets (many 4-point lines and concyclic quadruples), plus the (C1) facts
(blocks pairwise share <= 2 points, lines pairwise <= 1).  Uses count7.blocks_of (own code)."""
import itertools, random
from math import comb
from count7 import blocks_of
def check(pts):
    n = len(pts); bl = blocks_of(pts)
    circles = sum(1 for k in bl if k[0] == "C"); ell = sum(1 for k in bl if k[0] == "L")
    D = sum(comb(len(v), 3) - 1 for v in bl.values() if len(v) >= 4)
    assert comb(n, 3) - D - ell == circles, (pts, circles, D, ell)
    assert sum(comb(len(v), 3) for v in bl.values()) == comb(n, 3)
    vals = list(bl.values()); keys = list(bl.keys())
    for i, j in itertools.combinations(range(len(vals)), 2):
        assert len(vals[i] & vals[j]) <= 2
        if keys[i][0] == "L" and keys[j][0] == "L": assert len(vals[i] & vals[j]) <= 1
    return circles
grid = [(x, y) for x in range(4) for y in range(4)]
cnt = 0
for sub in itertools.combinations(grid, 7):
    check(sub); cnt += 1
print("4x4 grid: identity + (C1) verified on", cnt, "7-subsets")
random.seed(7); cnt = 0
for _ in range(3000):
    n = random.choice([6, 7, 8, 9])
    pool = [(x, y) for x in range(5) for y in range(5)]
    pts = random.sample(pool, n); check(pts); cnt += 1
print("random 5x5-grid subsets: identity + (C1) verified on", cnt, "sets")
