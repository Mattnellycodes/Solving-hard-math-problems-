"""Analyse the enumerated case-A classes: verify the forced shape (6-block, R-block, 18 (2,2)-blocks),
extract the six matchings M_{rr'} of the 6-block induced by the pairs of R, report whether complementary
pairs carry the same matching, and run all post-filters."""
import sys, json, itertools
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat')
from filters import full_report, bundle_violations
data = json.load(open(sys.argv[1]))
for i, rec in enumerate(data):
    F = [frozenset(B) for B in rec['blocks']]; L = [frozenset(S) for S in rec['lines']]
    six = [B for B in F if len(B) == 6]; assert len(six) == 1; S6 = six[0]; R = frozenset(range(10)) - S6
    n22 = [B for B in F if len(B) == 4 and len(B & S6) == 2]
    others = [B for B in F if B != S6 and B not in n22]
    M = {}
    for r, r2 in itertools.combinations(sorted(R), 2):
        M[(r, r2)] = frozenset(frozenset(B & S6) for B in n22 if {r, r2} <= B)
    comp = [((6,7),(8,9)), ((6,8),(7,9)), ((6,9),(7,8))]
    same = [M[a] == M[b] for a, b in comp]
    distinct = len(set(M.values()))
    print(f"class {i}: count={rec['count']} copies={rec['copies']} |S6|=6 R-block={R in F} n22={len(n22)} other blocks={[sorted(B) for B in others]} "
          f"lines={sorted(len(S) for S in L)} distinct matchings={distinct} complementary-equal={same}")
    rep = full_report(F, L)
    bv = bundle_violations(F, L)
    print(f"   bundle violations: {len(bv)}")
