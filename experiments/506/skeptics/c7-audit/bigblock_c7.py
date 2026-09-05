"""Skeptic c7-audit: check Lemmas A, B (big blocks) and the Lemma C inequality on n = 7 by brute force
over ALL labelled (C1)-families F of blocks of size 4..6 (no caps at all).
  Lemma A (block of size 6): circles >= f(7) = 13.
  Lemma B (largest block 5): circles >= 5^2 - 5 + 1 - 3*2 = 15.
  Lemma C: with B a largest block, |B| = m, r = 7 - m, N2 = #blocks (any size, incl. 3-blocks)
           meeting B in exactly 2 points:  N2 >= r*C(m,2) - m*r*(r-1)/4.
ell_max is computed by exact max-clique branch and bound over F u uncovered triples (pairwise <= 1).
"""
import itertools
from math import comb
n = 7; pts = range(n)
S = {k: [frozenset(c) for c in itertools.combinations(pts, k)] for k in (3, 4, 5, 6)}
cands = S[4] + S[5] + S[6]
fams = []
def dfs(start, F):
    fams.append(tuple(F))
    for j in range(start, len(cands)):
        if all(len(cands[j] & A) <= 2 for A in F):
            F.append(cands[j]); dfs(j + 1, F); F.pop()
dfs(0, [])
def all_blocks(F):
    return list(F) + [t for t in S[3] if not any(t <= B for B in F)]
def ell_max(F):
    lc = all_blocks(F); best = [0]
    def rec(i, L):
        if len(L) + (len(lc) - i) <= best[0]: return
        if i == len(lc): best[0] = max(best[0], len(L)); return
        if all(len(lc[i] & m) <= 1 for m in L):
            L.append(lc[i]); rec(i + 1, L); L.pop()
        rec(i + 1, L)
    rec(0, []); return best[0]
minc = {5: 99, 6: 99}; nfam = {5: 0, 6: 0}; lemmaC_ok = True; worst = 99
for F in fams:
    if not F: continue
    m = max(len(B) for B in F)
    D = sum(comb(len(B), 3) - 1 for B in F)
    if m in (5, 6):
        c = 35 - D - ell_max(F); minc[m] = min(minc[m], c); nfam[m] += 1
    B = max(F, key=len); r = n - m
    N2 = sum(1 for X in all_blocks(F) if len(X & B) == 2)
    bound = r * comb(m, 2) - m * r * (r - 1) / 4
    if N2 < bound: lemmaC_ok = False; print("Lemma C VIOLATED:", [sorted(b) for b in F], N2, bound)
    worst = min(worst, N2 - bound)
print("families with largest block 6:", nfam[6], " min circles =", minc[6], " (Lemma A claims >= 13)")
print("families with largest block 5:", nfam[5], " min circles =", minc[5], " (Lemma B claims >= 15)")
print("Lemma C inequality holds for all", len(fams) - 1, "non-empty families:", lemmaC_ok, "; min slack =", worst)
