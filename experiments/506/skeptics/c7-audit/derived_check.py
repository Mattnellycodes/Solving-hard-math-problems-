"""Skeptic c7-audit: purely combinatorial check of the reduction 'Fano-complement structure ->
complete quadrilateral + three diagonal quadruples at every point'.
Input: the unique surviving block family from enum_c7.py (7 four-blocks = complements of Fano lines).
For each point p we form the derived lines {B \\ {p} : p in B} and the 'circles' {B : p not in B}
and verify:
  (a) 4 derived lines, each of 3 points, pairwise meeting in exactly one point, every other point
      on exactly two of them  -> a complete quadrilateral with vertex P_ij = L_i ∩ L_j;
  (b) the three blocks avoiding p are exactly Q \\ {P_ij, P_kl} for the three opposite pairs,
      i.e. the three diagonal quadruples of the Angle Lemma.
Also re-derives, from scratch, that 7 four-subsets of a 7-set pairwise sharing <= 2 points must be
the Fano complements (complements are 7 triples pairwise sharing <= 1 point covering all 21 pairs).
"""
import itertools
F = [frozenset(b) for b in [[0,1,2,3],[0,1,4,5],[0,2,4,6],[0,3,5,6],[1,2,5,6],[1,3,4,6],[2,3,4,5]]]
pts = set(range(7))
comp = [pts - B for B in F]
pairs = set()
for T in comp:
    for pr in itertools.combinations(sorted(T), 2):
        assert pr not in pairs; pairs.add(pr)
print("complements are 7 triples covering each of the 21 pairs exactly once (STS(7) = Fano):", len(pairs) == 21)
ok_all = True
for p in range(7):
    Q = pts - {p}
    lines = [B - {p} for B in F if p in B]
    circ = [B for B in F if p not in B]
    a = (len(lines) == 4 and all(len(l) == 3 for l in lines)
         and all(len(l & m) == 1 for l, m in itertools.combinations(lines, 2))
         and all(sum(1 for l in lines if q in l) == 2 for q in Q))
    # vertex labelling
    V = {}
    for i, j in itertools.combinations(range(4), 2):
        (q,) = tuple(lines[i] & lines[j]); V[(i, j)] = q
    assert len(set(V.values())) == 6
    opp = [({V[(0,1)], V[(2,3)]}), ({V[(0,2)], V[(1,3)]}), ({V[(0,3)], V[(1,2)]})]
    diag = sorted(frozenset(Q - o) for o in opp)
    b = set(circ) == set(diag) and len(circ) == 3
    print(f"p={p}: quadrilateral={a}; blocks avoiding p = {[sorted(c) for c in circ]}; "
          f"opposite pairs = {[sorted(o) for o in opp]}; equal to the 3 diagonal quadruples: {b}")
    ok_all &= a and b
print("REDUCTION TO ANGLE LEMMA HOLDS AT EVERY POINT:", ok_all)
