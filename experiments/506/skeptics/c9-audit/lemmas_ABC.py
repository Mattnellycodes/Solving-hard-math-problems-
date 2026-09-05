"""Audit of Lemmas A, B, C (largest-block lemmas) for n = 9, and of the counting facts used.

Hand re-derivations are in the audit text; here:
 (1) numeric values of the three bounds for n = 9 under each line cap (ell <= 10 / 11 / 12);
 (2) empirical check of Lemma C's inequality  N2 >= r*C(m,2) - m*r(r-1)/4  (N2 = number of blocks
     meeting a largest block B in exactly 2 points) on random rational point sets (exact arithmetic),
     and of the resulting bound circles >= 1 + N2 - ell;
 (3) the counting step of the b5 = 0 case: b4 >= 16 forces a point of degree >= 8;
 (4) the elementary bounds ell <= 12 (pairs) and 'at most 8 four-blocks through a point'.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import comb, ceil
import random, sys
sys.path.insert(0, "/home/user/Solving-hard-math-problems-/experiments/506/skeptics/c9-audit")
from count9 import blocks

def lemmaA(n):            # block of size n-1: circles >= f(n)
    return comb(n - 1, 2) + 1 - (n - 1) // 2
def lemmaB(n, ell_cap):   # block of size m = n-2: |B(P)| = 1 + 2C(m,2) + m - 3t, t <= floor(m/2); circles = |B(P)| - ell
    m = n - 2
    return 1 + 2 * comb(m, 2) + m - 3 * (m // 2) - min(ell_cap, m)   # ell <= m shown in Lemma B; also ell<=cap
def lemmaB_weak(n, ell_cap):
    m = n - 2
    return 1 + 2 * comb(m, 2) + m - 3 * (m // 2) - ell_cap            # not even using ell <= m
def lemmaC(n, m, ell_cap):
    r = n - m
    return 1 + ceil(r * comb(m, 2) - m * r * (r - 1) / 4) - ell_cap

print("(1) n = 9 largest-block bounds on circles (need > 24 to exclude):")
for ell_cap, name in ((10, "strong o(9)=6"), (11, "SG only"), (12, "no cap")):
    print(f"   ell <= {ell_cap} ({name}): m=8 (Lemma A) >= {lemmaA(9)} ; m=7 (Lemma B) >= {lemmaB(9, ell_cap)} "
          f"(without ell<=m: {lemmaB_weak(9, ell_cap)}) ; m=6 (Lemma C) >= {lemmaC(9, 6, ell_cap)}")

print("(2) empirical check of Lemma C on random rational point sets (exact arithmetic):")
random.seed(1)
viol = 0; tested = 0; worst = None
for trial in range(400):
    n = random.choice([7, 8, 9, 10])
    # bias towards rich structures: points from a small grid or from a circle + grid
    P = set()
    while len(P) < n:
        if random.random() < 0.5:
            P.add((Fr(random.randint(-2, 2)), Fr(random.randint(-2, 2))))
        else:
            t = Fr(random.randint(-3, 3), random.randint(1, 3))
            P.add((Fr(1 - t * t, 1 + t * t), Fr(2 * t, 1 + t * t)))
    P = list(P)
    B = blocks(P)
    if any(len(s) == n for s in B.values()):
        continue
    tested += 1
    big = max(B.values(), key=len)
    m = len(big); r = n - m
    if m < 4: continue
    N2 = sum(1 for s in B.values() if len(s & big) == 2)
    bound = r * comb(m, 2) - Fr(m * r * (r - 1), 4)
    circles = sum(1 for k in B if k[0] != 0); ell = sum(1 for k in B if k[0] == 0)
    if N2 < bound or circles < 1 + N2 - ell:
        viol += 1; print("   VIOLATION", P, m, N2, bound, circles, ell)
    slack = N2 - bound
    if worst is None or slack < worst[0]: worst = (slack, n, m, N2, bound)
print(f"   tested {tested} non-degenerate sets, violations {viol}, smallest slack N2 - bound = {worst}")

print("(3) b5 = 0 counting: b4 >= 16 four-blocks on 9 points => sum of degrees >= 64 > 9*7 => some point in >= 8 four-blocks:",
      4 * 16 > 9 * 7)
print("    with D + ell >= 60 and ell <= 12: D = 3*b4 >= 48 => b4 >= 16:", 3 * 15 < 48)
print("(4) ell <= floor(36/3) = 12 lines on 9 points (pairwise <= 1 common point, each covers >= 3 pairs):", 36 // 3)
print("    four-blocks through a point: derived triples on 8 points pairwise sharing <= 1 point: each point in <= 3 => <= 8")
