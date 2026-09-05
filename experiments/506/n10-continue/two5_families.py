"""Family-level check for every F-class of the two-disjoint-5-block skeleton: real kernel dimension of the
angle system, torsion exponent, and the number of torsion labellings with distinct points on each circle."""
import sys, json, itertools
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-continue')
import sympy as sp
from analyse_two5_own import torsion_solutions
d = json.load(open('runs/enum_7.json'))
summary = {}
for k, rec in enumerate(d['F_classes']):
    F = [frozenset(B) for B in rec['blocks']]
    fives = sorted([B for B in F if len(B) == 5], key=sorted); S, R = sorted(fives[0]), sorted(fives[1])
    idx = {p: i for i, p in enumerate(S + R)}
    fours = [B for B in F if len(B) == 4]
    M = [[0] * 10 for _ in fours]
    for j, B in enumerate(fours):
        for p in B: M[j][idx[p]] = 1 if p in S else -1
    Ms = sp.Matrix(M); kd = len(Ms.nullspace())
    if kd != 1:
        print(f"F-class {k}: b4={rec['b4']} kernel dim {kd} -> NOT DECIDED by torsion method"); summary[k] = ('kernel>1', kd); continue
    sols, e, r = torsion_solutions(M)
    norm = {tuple((x - s[0]) % e for x in s) for s in sols}
    good = sorted(s for s in norm if len(set(s[:5])) == 5 and len(set(s[5:])) == 5)
    print(f"F-class {k}: b4={rec['b4']} line sets={rec['labelled_line_sets']} e={e} torsion={len(sols)} distinct-point labellings={len(good)}" + (f"  e.g. {good[0]}" if good else "  => 4-BLOCK FAMILY NOT REALISABLE (case i); cases ii/iii excluded by kernel"))
    summary[k] = ('labellings', len(good))
print("F-classes with at least one distinct-point torsion labelling:", [k for k, v in summary.items() if v[0] == 'labellings' and v[1] > 0])
print("F-classes with kernel dim > 1:", [k for k, v in summary.items() if v[0] == 'kernel>1'])
