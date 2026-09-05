"""Hostile-referee certificate for the Angle Lemma (REPORT.md Lemma 2.3) -- light version.

Lemma. Four lines L1..L4 in R^2, pairwise non-parallel, no three concurrent, P_ij = L_i ∩ L_j.
The three quadruples  Q1 = {P12,P34,P13,P24},  Q2 = {P12,P34,P14,P23},  Q3 = {P13,P24,P14,P23}
are never all concyclic.

Strategy (mine, independent of the theory agent's script):
 (S) The set {Q1,Q2,Q3} is invariant under relabelling the lines (each Q_k is the complement of a pair of
     'diagonal' vertices {P_ab, P_cd}, {ab|cd} a split of {1,2,3,4}).  Checked below by brute force.
 (1) Chart 1 -- no line perpendicular to L1.  Similarity normalisation: L1: y = 0, P12 = origin,
     L3 has y-intercept 1 (L3 does not pass through P12).   L2: y = a x, L3: y = b x + 1, L4: y = c x + d,
     a, b, c finite (no line perpendicular to L1), pairwise distinct, all non-zero;  d != 0 (L1,L2,L4 not
     concurrent);  c != b d (L1,L3,L4 not concurrent);  P23 not on L4 (L2,L3,L4 not concurrent).
     We factor the three concyclicity determinants.  Every factor other than one 'essential' factor
     A_k must be a non-degeneracy polynomial (so it cannot vanish); the essential factors satisfy
     A1 - A2 = 2(b - c)  (or similar) => any two of the three conditions force two lines parallel.
 (2) Chart 2 -- some line is perpendicular to L1; by (S) WLOG it is L2: L2: x = 0 (then L3, L4 are not
     vertical, else parallel to L2).  L3: y = b x + 1, L4: y = c x + d, b, c non-zero distinct, d != 0,
     c != b d, P23 not on L4.  Factor again; essential factors should be  bc + 1  (Q1, Q2) and  bc - 1  (Q3),
     which cannot vanish simultaneously.
 (3) Tiny Groebner check that the essential factors + non-parallelism generate the unit ideal.
 (4) Numerical sanity: random non-degenerate quadrilaterals; and in chart 2 the exact 2-condition family
     (b c = -1) with Q3 evaluated on it.
Implementation note: points are kept as (X, Y, W) with polynomial X, Y, W (P = (X/W, Y/W)); the row
[x, y, x^2+y^2, 1] is replaced by the polynomial row [X W, Y W, X^2+Y^2, W^2] = W^2 * original row, so
the determinant is a polynomial equal to (prod W^2) * (concyclicity determinant).  W factors are
non-degeneracy factors (denominators of intersection points = non-parallelism), so nothing is lost.
"""
import itertools, random
import sympy as sp

a, b, c, d, t = sp.symbols('a b c d t')
SPLITS = {1: ((1, 4), (2, 3)), 2: ((1, 3), (2, 4)), 3: ((1, 2), (3, 4))}   # Q_k = complement of these diagonal pairs


def symmetry_check():
    pairs = [frozenset(p) for p in itertools.combinations(range(1, 5), 2)]
    quads = set()
    for k, (p, q) in SPLITS.items():
        quads.add(frozenset(pairs) - {frozenset(p), frozenset(q)})
    assert len(quads) == 3
    ok = True
    for s in itertools.permutations(range(1, 5)):
        m = dict(zip(range(1, 5), s))
        img = {frozenset(frozenset(m[i] for i in pr) for pr in Q) for Q in quads}
        ok &= (img == quads)
    return ok


def inter(L1, L2):
    """lines as (A, B, C) meaning A x + B y = C; returns homogeneous (X, Y, W) with point (X/W, Y/W)."""
    (A1, B1, C1), (A2, B2, C2) = L1, L2
    W = sp.expand(A1 * B2 - A2 * B1)
    X = sp.expand(C1 * B2 - C2 * B1)
    Y = sp.expand(A1 * C2 - A2 * C1)
    return (X, Y, W)


def concyc_poly(pts):
    rows = [[X * W, Y * W, X * X + Y * Y, W * W] for (X, Y, W) in pts]
    return sp.expand(sp.Matrix(rows).det(method='berkowitz'))


def conditions(Lines):
    P = {(i, j): inter(Lines[i], Lines[j]) for i, j in itertools.combinations(range(1, 5), 2)}
    Qs = {}
    for k, (pr, qr) in SPLITS.items():
        quad = [P[ij] for ij in itertools.combinations(range(1, 5), 2) if ij not in (pr, qr)]
        Qs[k] = sp.factor(concyc_poly(quad))
    return P, Qs


def analyse(chart_name, Lines, nondeg_factors):
    print(f"=== {chart_name}", flush=True)
    P, Qs = conditions(Lines)
    for ij, (X, Y, W) in P.items():
        print(f"   P{ij[0]}{ij[1]} = ({X}, {Y}) / ({W})")
    nd_set = {sp.factor(f) for f in nondeg_factors}
    essential = {}
    for k, Q in Qs.items():
        fac = sp.factor_list(Q)
        others = []
        for f, e in fac[1]:
            f = sp.factor(f)
            if f in nd_set or sp.factor(-f) in nd_set:
                continue
            others.append((f, e))
        print(f"   Q{k} (times denominators) = {Q}")
        assert len(others) == 1, f"unexpected factor structure: {others}"
        essential[k] = sp.expand(others[0][0])
        print(f"      essential factor A{k} = {essential[k]}   (all other factors are non-degeneracy polynomials)")
    return essential


def main():
    print("S4 symmetry of {Q1,Q2,Q3} under relabelling lines:", symmetry_check(), flush=True)
    # chart 1: y = m x + k  ->  -m x + y = k
    L = {1: (sp.Integer(0), sp.Integer(1), sp.Integer(0)), 2: (-a, sp.Integer(1), sp.Integer(0)),
         3: (-b, sp.Integer(1), sp.Integer(1)), 4: (-c, sp.Integer(1), d)}
    P = {(i, j): inter(L[i], L[j]) for i, j in itertools.combinations(range(1, 5), 2)}
    X23, Y23, W23 = P[(2, 3)]
    conc234 = sp.expand(Y23 - c * X23 - d * W23)          # P23 on L4  <=>  y = c x + d  (times W23)
    nondeg1 = [a, b, c, a - b, a - c, b - c, d, c - b * d, conc234]
    print("chart-1 non-degeneracy polynomials:", [sp.factor(f) for f in nondeg1], flush=True)
    E1 = analyse("chart 1: L1: y=0, L2: y=ax, L3: y=bx+1, L4: y=cx+d", L, nondeg1)
    print("   A1 - A2 =", sp.factor(E1[1] - E1[2]), ";  A1 - A3 =", sp.factor(E1[1] - E1[3]),
          ";  A2 - A3 =", sp.factor(E1[2] - E1[3]), flush=True)
    G = sp.groebner([E1[1], E1[2], E1[3], t * (a - b) * (a - c) * (b - c) - 1], a, b, c, d, t, order='grevlex')
    print("   Groebner(A1, A2, A3, non-parallel) == {1}:", list(G) == [1], flush=True)
    for (i, j) in [(1, 2), (1, 3), (2, 3)]:
        Gij = sp.groebner([E1[i], E1[j], t * (a - b) * (a - c) * (b - c) - 1], a, b, c, d, t, order='grevlex')
        print(f"   Groebner(A{i}, A{j}, non-parallel) == {{1}} (any TWO conditions impossible in chart 1):", list(Gij) == [1], flush=True)
    # chart 2: L2 vertical x = 0
    L2 = {1: (sp.Integer(0), sp.Integer(1), sp.Integer(0)), 2: (sp.Integer(1), sp.Integer(0), sp.Integer(0)),
          3: (-b, sp.Integer(1), sp.Integer(1)), 4: (-c, sp.Integer(1), d)}
    P2 = {(i, j): inter(L2[i], L2[j]) for i, j in itertools.combinations(range(1, 5), 2)}
    X23, Y23, W23 = P2[(2, 3)]
    conc234b = sp.expand(Y23 - c * X23 - d * W23)
    nondeg2 = [b, c, b - c, d, c - b * d, conc234b]
    print("chart-2 non-degeneracy polynomials:", [sp.factor(f) for f in nondeg2], flush=True)
    E2 = analyse("chart 2: L1: y=0, L2: x=0, L3: y=bx+1, L4: y=cx+d", L2, nondeg2)
    G2 = sp.groebner([E2[1], E2[2], E2[3]], b, c, d, order='grevlex')
    print("   Groebner(A1, A2, A3) in chart 2:", list(G2), " == {1}:", list(G2) == [1], flush=True)
    G2_12 = sp.groebner([E2[1], E2[2]], b, c, d, order='grevlex')
    print("   Groebner(A1, A2) in chart 2 (the two-condition family):", list(G2_12), flush=True)
    # numeric sanity
    random.seed(3)
    Q1f = sp.lambdify((a, b, c, d), [E1[1], E1[2], E1[3]], 'math')
    nd1 = sp.lambdify((a, b, c, d), sp.prod(nondeg1), 'math')
    best = None
    for _ in range(50000):
        v = [random.uniform(-4, 4) for _ in range(4)]
        if abs(nd1(*v)) < 1e-4:
            continue
        m = max(abs(x) for x in Q1f(*v))
        if best is None or m < best[0]:
            best = (m, v)
    print("   numeric chart 1: min over 50000 random non-degenerate quadrilaterals of max_k |A_k| =", best)
    Q2f = sp.lambdify((b, c, d), [E2[1], E2[2], E2[3]], 'math')
    worst = None
    for _ in range(10000):
        bb = random.uniform(-4, 4)
        if abs(bb) < 1e-3:
            continue
        cc = -1 / bb; dd = random.uniform(-4, 4)
        vals = Q2f(bb, cc, dd)
        if worst is None or abs(vals[2]) < worst[0]:
            worst = (abs(vals[2]), (bb, cc, dd), vals)
    print("   numeric chart 2 on the exact family L3 ⊥ L4 (A1 = A2 = 0): min |A3| =", worst)


if __name__ == "__main__":
    main()
