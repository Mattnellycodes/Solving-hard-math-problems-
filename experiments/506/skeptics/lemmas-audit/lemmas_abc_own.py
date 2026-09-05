"""lemmas-audit: brute-force check of Lemmas A, B, C on ALL abstract block structures of small n.
An abstract structure = family F of blocks of size >= 4 pairwise sharing <= 2 points (every triple not
inside a member of F is a 3-block).  Lemmas A, B, C are purely combinatorial statements about such
structures plus the line rules (lines pairwise share <= 1 point; a line of size >= 4 is in F; a line of
size 3 is an uncovered triple).  We compute for every structure (from enum_own with capmode 0 and
target = C(n,3), i.e. ALL families) the largest block B, m = |B|, r = n - m,
   N2 = number of blocks (any size, 3-blocks included) meeting B in exactly 2 points,
   ell_max (no caps, pure line rules), circles_min = C(n,3) - D - ell_max,
and test:
   Lemma C: N2 >= r*C(m,2) - m*r*(r-1)/4         (for every structure, r >= 1)
   Lemma B: if m = n-2: circles_min >= m^2 - m + 1 - 3*floor(m/2)   and ell_max <= m
   Lemma A: if m = n-1: circles_min >= f(n)     and equality only if x lies on floor((n-1)/2) lines
usage: python3 lemmas_abc_own.py n file.out
"""
import sys, itertools
from math import comb

def f(n): return comb(n - 1, 2) + 1 - (n - 1) // 2

def decode(mask, n): return frozenset(i for i in range(n) if mask >> i & 1)

def all_blocks(F, n):
    """blocks of every size: F plus the uncovered triples."""
    return list(F) + [frozenset(t) for t in itertools.combinations(range(n), 3) if not any(frozenset(t) <= b for b in F)]

def ell_max(F, n, capl):
    cands = all_blocks(F, n)
    best = [0]
    def rec(start, chosen, pairs):
        best[0] = max(best[0], len(chosen))
        avail = [j for j in range(start, len(cands)) if all(len(cands[j] & c) <= 1 for c in chosen) and pairs + comb(len(cands[j]), 2) <= capl]
        if len(chosen) + len(avail) <= best[0]:
            return
        for j in avail:
            rec(j + 1, chosen + [cands[j]], pairs + comb(len(cands[j]), 2))
    rec(0, [], 0)
    return best[0]

def main():
    n = int(sys.argv[1]); fn = sys.argv[2]
    N3 = comb(n, 3); capl = comb(n, 2)
    stats = {"A": [0, 0], "B": [0, 0], "C": [0, 0]}
    worstC = None; tightC = 0
    seen = 0
    for line in open(fn):
        if not line.startswith("REC"): continue
        head, masks = line.split(":")
        F = [decode(int(m), n) for m in masks.split()]
        if not F:
            continue
        seen += 1
        D = sum(comb(len(b), 3) - 1 for b in F)
        B = max(F, key=len); m = len(B); r = n - m
        blocks = all_blocks(F, n)
        assert sum(comb(len(b), 3) for b in blocks) == N3
        N2 = sum(1 for b in blocks if len(b & B) == 2)
        boundC = r * comb(m, 2) - m * r * (r - 1) / 4
        stats["C"][0] += 1
        if N2 < boundC:
            stats["C"][1] += 1; print("LEMMA C VIOLATED", [sorted(b) for b in F], N2, boundC)
        if N2 == boundC: tightC += 1
        if worstC is None or N2 - boundC < worstC[0]:
            worstC = (N2 - boundC, [sorted(b) for b in F])
        if r <= 2:
            e = ell_max(F, n, capl)
            circles = N3 - D - e
            if r == 2:
                stats["B"][0] += 1
                if circles < m * m - m + 1 - 3 * (m // 2) or e > m:
                    stats["B"][1] += 1; print("LEMMA B VIOLATED", [sorted(b) for b in F], circles, e)
            if r == 1:
                stats["A"][0] += 1
                if circles < f(n):
                    stats["A"][1] += 1; print("LEMMA A VIOLATED", [sorted(b) for b in F], circles)
    print(f"n={n}: structures checked={seen}; Lemma C: {stats['C'][0]} checked, {stats['C'][1]} violations, "
          f"tight in {tightC}, min slack {worstC[0]} at {worstC[1]}")
    print(f"      Lemma B (m=n-2): {stats['B'][0]} checked, {stats['B'][1]} violations; Lemma A (m=n-1): {stats['A'][0]} checked, {stats['A'][1]} violations")

if __name__ == "__main__":
    main()
