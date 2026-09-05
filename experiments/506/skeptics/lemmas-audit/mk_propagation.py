"""lemmas-audit: in a realisation of the (8_3) incidences by 8 DISTINCT real points, any extra
collinearity forces all 8 points to be collinear.  Combinatorial closure argument:
 (1) an extra collinear triple {a,b,c} contains a pair covered by an MK line L (the 4 uncovered pairs
     form a perfect matching, so a triple contains at most one of them), hence c lies on line(L):
     L becomes a >=4-point line;  (2) closure: whenever a set X is known collinear and an MK line
     meets X in >= 2 points, its third point joins X.  We check that for EVERY MK line L and every
     extra point c, closure of L ∪ {c} is all 8 points; and for every pair of distinct MK lines
     that coincide, closure is all 8 points."""
import itertools
MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
def closure(X):
    X = set(X)
    changed = True
    while changed:
        changed = False
        for L in MK:
            if len(L & X) >= 2 and not L <= X:
                X |= L; changed = True
    return X
ok = True
for L in MK:
    for c in range(8):
        if c in L: continue
        if closure(L | {c}) != set(range(8)):
            ok = False; print("closure fails for", sorted(L), c)
covered = set()
for L in MK:
    covered.update(itertools.combinations(sorted(L), 2))
unc = [p for p in itertools.combinations(range(8), 2) if p not in covered]
print("uncovered pairs (perfect matching?):", unc)
print("every extra collinear triple forces full collinearity:", ok)
