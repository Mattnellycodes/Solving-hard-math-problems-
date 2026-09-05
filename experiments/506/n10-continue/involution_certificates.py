"""Explicit certificates for the single-6-block F-classes: lex Groebner basis of the saturated involution ideal."""
import sys, json, itertools
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-continue')
import sympy as sp
from analyse_six_own import t3, t4, t5, w
d = json.load(open('runs/enum_1.json'))
for k, rec in enumerate(d['F_classes']):
    F = [frozenset(B) for B in rec['blocks']]
    A = sorted(next(B for B in F if len(B) == 6)); R = [p for p in range(10) if p not in A]
    lab = {A[0]: (sp.Integer(0), sp.Integer(1)), A[1]: (sp.Integer(1), sp.Integer(1)), A[2]: (sp.Integer(1), sp.Integer(0)),
           A[3]: (t3, sp.Integer(1)), A[4]: (t4, sp.Integer(1)), A[5]: (t5, sp.Integer(1))}
    conds = []
    for r, r2 in itertools.combinations(R, 2):
        pairs = [tuple(sorted(B & set(A))) for B in F if {r, r2} <= B and len(B) == 4 and len(B & set(A)) == 2]
        assert len(pairs) == 3 and len(set(itertools.chain.from_iterable(pairs))) == 6
        rows = [[lab[a][0] * lab[b][0], lab[a][0] * lab[b][1] + lab[b][0] * lab[a][1], lab[a][1] * lab[b][1]] for a, b in pairs]
        conds.append(sp.expand(sp.Matrix(rows).det()))
    dist = sp.Integer(1)
    for a, b in itertools.combinations(A, 2):
        dist *= (lab[a][0] * lab[b][1] - lab[b][0] * lab[a][1])
    Gl = sp.groebner(conds + [sp.expand(dist * w - 1)], w, t3, t4, t5, order='lex')
    sat = [g for g in Gl.exprs if not g.has(w)]
    print(f"F-class {k} (line sets: {rec['labelled_line_sets']}): involution conditions {conds}")
    print(f"    saturated ideal, lex basis (t3 > t4 > t5): {sat}")
