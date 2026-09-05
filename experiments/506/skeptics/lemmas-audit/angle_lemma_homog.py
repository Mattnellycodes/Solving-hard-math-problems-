"""lemmas-audit: fraction-free algebraic certificate of the Angle Lemma in the two charts that
cover all cases (see angle_lemma_light.py for the reduction):
  chart A: L1: y = 0, L2: y = m2 x (through P12 = origin), L3: y = m3 x + c3, L4: y = m4 x + c4
           (no line perpendicular to L1; this is the theory agent's chart);
  chart B: L1: y = 0, L2: x = 0 (perpendicular), L3: y = m3 x + c3, L4: y = m4 x + c4
           (the chart the theory agent's certificate does not cover).
Lines as vectors (a, b, c) meaning a x + b y + c = 0; vertex P_ij = L_i x L_j (homogeneous).
Four homogeneous points (X:Y:Z) are concyclic iff det[[X Z, Y Z, X^2 + Y^2, Z^2]] = 0 (rows).
Non-degeneracy: pairwise non-parallel (a_i b_j - a_j b_i != 0) and no three concurrent
(det of the three line vectors != 0); this makes the six vertices distinct.
Method: factor each concyclicity polynomial; strip non-degeneracy factors; show the essential
factors are jointly inconsistent with non-degeneracy by a small Groebner computation.
"""
import itertools, sympy as sp

def cross(u, v):
    return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])

def concyc(pts):
    rows = [[X * Z, Y * Z, X * X + Y * Y, Z * Z] for (X, Y, Z) in pts]
    return sp.expand(sp.Matrix(rows).det())

def run(name, lines, syms):
    P = {(i, j): cross(lines[i], lines[j]) for i, j in itertools.combinations(range(1, 5), 2)}
    nd = []
    for i, j in itertools.combinations(range(1, 5), 2):
        nd.append(sp.expand(lines[i][0] * lines[j][1] - lines[j][0] * lines[i][1]))
    for i, j, k in itertools.combinations(range(1, 5), 3):
        nd.append(sp.expand(sp.Matrix([lines[i], lines[j], lines[k]]).det()))
    ndf = set()
    for q in nd:
        if q == 0:
            print(f"  chart {name}: a non-degeneracy condition is identically 0 -> chart invalid"); return None
        for fac, _ in sp.factor_list(q)[1]:
            ndf.add(sp.expand(fac)); ndf.add(sp.expand(-fac))
    Q = [concyc([P[(1, 2)], P[(3, 4)], P[(1, 3)], P[(2, 4)]]),
         concyc([P[(1, 2)], P[(3, 4)], P[(1, 4)], P[(2, 3)]]),
         concyc([P[(1, 3)], P[(2, 4)], P[(1, 4)], P[(2, 3)]])]
    print(f"--- chart {name}; non-degeneracy factors: {sorted(str(f) for f in ndf if str(f)[0] != '-')}")
    ess = []
    for i, q in enumerate(Q, 1):
        fl = sp.factor_list(q)
        e = []; d = []
        for fac, mult in fl[1]:
            (d if sp.expand(fac) in ndf else e).append(fac)
        print(f"  Q{i}: degenerate factors {d}; essential {e}")
        ess.append(e)
    E = [sp.expand(x) for e in ess for x in e]
    G = sp.groebner(E, *syms, order="lex")
    print("  lex Groebner basis of essential factors:", list(G))
    t = sp.Symbol("t")
    prod = sp.Integer(1)
    for fct in {f for f in ndf if str(f)[0] != '-'}:
        prod *= fct
    G2 = sp.groebner(E + [sp.expand(t * prod - 1)], *syms, t, order="grevlex")
    ok = list(G2) == [1]
    print("  saturated by non-degeneracy -> unit ideal:", ok)
    return ok

if __name__ == "__main__":
    m2, m3, m4, c3, c4 = sp.symbols("m2 m3 m4 c3 c4")
    A = {1: (0, 1, 0), 2: (-m2, 1, 0), 3: (-m3, 1, -c3), 4: (-m4, 1, -c4)}
    okA = run("A (no line perpendicular to L1)", A, [m2, m3, m4, c3, c4])
    B = {1: (0, 1, 0), 2: (1, 0, 0), 3: (-m3, 1, -c3), 4: (-m4, 1, -c4)}
    okB = run("B (L1 perpendicular to L2)", B, [m3, m4, c3, c4])
    print("ANGLE LEMMA CERTIFIED IN BOTH CHARTS:", bool(okA and okB))
    # Illustrate why chart B matters: in chart B, Q1 and Q2 essential factors alone
    # have real solutions (L1 ⊥ L2 with suitable L3, L4) -- the third quadruple is needed.
