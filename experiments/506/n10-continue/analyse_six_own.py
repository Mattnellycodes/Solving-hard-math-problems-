"""Exact non-realisability of 4-block families containing a 6-block A (own derivation and code).

Lemma (pencil involution).  Sphere model: the circle A lies in a plane π; r, r' ∈ P off A.  Every circle through
r, r' is a plane section through the line rr', which meets π in one point X (possibly at infinity of π; X ∉ A,
because otherwise every block through r, r' would contain X, impossible for three blocks cutting A in disjoint
pairs).  So the pairs {a, b} ⊂ A cut out by blocks through r, r' are pairs of the projective involution of the
conic A with centre X (chords through X).  If three blocks {r, r', a_i, b_i} pair up all six points of A, the
three pairs belong to one involution:  with A ≅ P^1 (parameters (x_i : y_i)) an involution is a symmetric form
α x x' + β (x y' + x' y) + γ y y' = 0, so   det[[x_a x_b, x_a y_b + x_b y_a, y_a y_b]]_{3 pairs} = 0.
PGL(2, R) is sharply 3-transitive: put A_0 = (0:1), A_1 = (1:1), A_2 = (1:0); the others are (t_i : 1) with t_i
finite (distinct from A_2) and pairwise distinct, ≠ 0, 1.  Groebner basis (1) of the determinant conditions
saturated by the distinctness product  =>  no realisation over C, hence none over R, whatever the lines.
Usage: python3 analyse_six_own.py enum_file.json"""
import sys, json, itertools, time
import sympy as sp
t3, t4, t5, w = sp.symbols('t3 t4 t5 w')
def analyse(F):
    F = [frozenset(B) for B in F]
    A = sorted(next(B for B in F if len(B) == 6)); R = [p for p in range(10) if p not in A]
    lab = {A[0]: (sp.Integer(0), sp.Integer(1)), A[1]: (sp.Integer(1), sp.Integer(1)), A[2]: (sp.Integer(1), sp.Integer(0)),
           A[3]: (t3, sp.Integer(1)), A[4]: (t4, sp.Integer(1)), A[5]: (t5, sp.Integer(1))}
    conds = []; info = []
    for r, r2 in itertools.combinations(R, 2):
        pairs = [tuple(sorted(B & set(A))) for B in F if {r, r2} <= B and len(B) == 4 and len(B & set(A)) == 2]
        if len(pairs) == 3 and len(set(itertools.chain.from_iterable(pairs))) == 6:
            rows = [[lab[a][0] * lab[b][0], lab[a][0] * lab[b][1] + lab[b][0] * lab[a][1], lab[a][1] * lab[b][1]] for a, b in pairs]
            conds.append(sp.expand(sp.Matrix(rows).det())); info.append(((r, r2), pairs))
    dist = sp.Integer(1)
    for a, b in itertools.combinations(A, 2):
        dist *= (lab[a][0] * lab[b][1] - lab[b][0] * lab[a][1])
    t0 = time.time()
    G = sp.groebner(conds + [sp.expand(dist * w - 1)], w, t3, t4, t5, order='grevlex')
    unit = list(G.exprs) == [1]
    real_sols = None
    if not unit:
        # eliminate w (lex, w first): generators of the saturated ideal in Q[t3, t4, t5]; solve exactly
        Gl = sp.groebner(conds + [sp.expand(dist * w - 1)], w, t3, t4, t5, order='lex')
        sat = [g for g in Gl.exprs if not g.has(w)]
        sols = sp.solve(sat, [t3, t4, t5], dict=True)
        real_sols = [s for s in sols if all(sp.simplify(sp.im(val)) == 0 for val in s.values()) and len(s) == 3]
        param = [s for s in sols if len(s) < 3]
        real_sols = (real_sols, param, sols)
    return unit, info, conds, time.time() - t0, real_sols
if __name__ == '__main__':
    d = json.load(open(sys.argv[1]))
    recs = d['F_classes'] if 'F_classes' in d else d
    for k, rec in enumerate(recs):
        F = rec['blocks']
        if not any(len(B) == 6 for B in F): continue
        unit, info, conds, dt, rs = analyse(F)
        print(f"F-class {k}: 6-block matchings per R-pair: {[(rr, m) for rr, m in info]}")
        if unit:
            verdict = "4-BLOCK FAMILY NOT REALISABLE over C (any lines)"
        else:
            real_sols, param, sols = rs
            if not param and not real_sols:
                verdict = f"saturated ideal zero-dimensional with {len(sols)} complex solutions, none real => 4-BLOCK FAMILY NOT REALISABLE over R (any lines)"
            else:
                verdict = f"REAL SOLUTIONS MAY EXIST: isolated real {real_sols}, parametric {param} => needs the lines analysis"
        print(f"   {len(conds)} involution conditions; saturated Groebner basis = (1): {unit}  [{dt:.1f}s]")
        if not unit: print("   all solutions of the saturated system:", rs[2])
        print("   =>", verdict)
