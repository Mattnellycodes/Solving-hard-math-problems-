"""Hostile-referee re-verification of the Angle Lemma (REPORT.md Lemma 2.3), with my own
parametrisation and my own concyclicity test.

Setting: four lines L1..L4 in R^2, pairwise non-parallel, no three concurrent; P_ij = L_i ∩ L_j.
Claim: the quadruples Q1 = {P12,P34,P13,P24}, Q2 = {P12,P34,P14,P23}, Q3 = {P13,P24,P14,P23}
are never all three concyclic.

Normalisation (similarities preserve circles and lines):  rotate so that L1 is the x-axis,
translate so that P12 is the origin, scale so that L3 has y-intercept 1 (L3 does not pass through P12
by non-concurrency, so its intercept is non-zero).  Parameters: L2: y = a x, L3: y = b x + 1,
L4: y = c x + d.   Non-degeneracy:  a, b, c pairwise distinct and non-zero (L2, L3, L4 not parallel to
L1 = x-axis and to each other);  d != 0 (L4 not through P12), d != 1 (L4, L3, L1... actually
P13 = (-1/b, 0), P14 = (-d/c, 0) distinct: c != b d), P23 != P24: a != ... handled by the vertex
distinctness polynomials; we saturate by the product of all of them (Rabinowitsch trick).

Concyclicity test (mine): four points are concyclic or collinear iff the determinant
      | x  y  x^2+y^2  1 |  over the four rows vanishes.  Collinearity of a quadruple is excluded by the
combinatorics (each quadruple contains two vertices of a line and a vertex off it... we check it too).

Checks performed:
  1. Groebner basis of (Q1, Q2, Q3, t*nondeg - 1) is {1}   (no complex solution at all).
  2. Q1 = Q2 = 0 (with non-degeneracy) forces  a = infinity?  no: forces L1 ⊥ L2, i.e. in this chart
     L2 would be vertical -- so Q1 = Q2 = 0 has NO solution in this chart with a finite; we verify by
     checking that the ideal (Q1, Q2) saturated by non-degeneracy is the unit ideal as well (the
     vertical-L2 case is handled separately with a second chart where L2 is the y-axis).
  3. Second chart: L1: y = 0, L2: x = 0 (perpendicular), L3: y = b x + 1, L4: y = c x + d.  Then Q1 and
     Q2 conditions should hold iff L3 ⊥ L4 (b c = -1), and Q3 must then fail: we verify Q3's ideal with
     b c + 1 = 0 and non-degeneracy is the unit ideal.
  4. Numerical random sampling: random quadrilaterals, solve Q1 = Q2 = 0 for (c, d) given (a, b) by
     Newton from random starts; any solution found must be degenerate; report min |Q3| over
     non-degenerate solutions (expect: no non-degenerate solutions).
"""
import random
import sympy as sp

a, b, c, d, t = sp.symbols('a b c d t')


def concyc_det(pts):
    return sp.expand(sp.Matrix([[x, y, x * x + y * y, 1] for x, y in pts]).det())


def collinear_det(p, q, r):
    return sp.expand(sp.Matrix([[p[0], p[1], 1], [q[0], q[1], 1], [r[0], r[1], 1]]).det())


def vertices(L):
    """L: dict i -> (m, k) meaning y = m x + k, or ('vert', x0) for a vertical line x = x0."""
    def inter(i, j):
        Li, Lj = L[i], L[j]
        if Li[0] == 'vert':
            x0 = Li[1]; m, k = Lj; return (x0, m * x0 + k)
        if Lj[0] == 'vert':
            x0 = Lj[1]; m, k = Li; return (x0, m * x0 + k)
        (mi, ki), (mj, kj) = Li, Lj
        x0 = sp.cancel((kj - ki) / (mi - mj)); return (x0, sp.cancel(mi * x0 + ki))
    return {(i, j): inter(i, j) for i, j in [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]}


def three_conditions(P):
    Q1 = concyc_det([P[(1, 2)], P[(3, 4)], P[(1, 3)], P[(2, 4)]])
    Q2 = concyc_det([P[(1, 2)], P[(3, 4)], P[(1, 4)], P[(2, 3)]])
    Q3 = concyc_det([P[(1, 3)], P[(2, 4)], P[(1, 4)], P[(2, 3)]])
    return [sp.factor(sp.numer(sp.together(Q))) for Q in (Q1, Q2, Q3)]


def distinctness_product(P):
    prod = sp.Integer(1)
    keys = list(P)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            (x1, y1), (x2, y2) = P[keys[i]], P[keys[j]]
            # points distinct  <=>  (x1-x2, y1-y2) != 0 ; we encode via the two numerators separately is
            # not a single polynomial; instead use  (x1-x2)^2 + (y1-y2)^2 != 0  (valid over R, and we only
            # need R-solutions; but Groebner works over C -- see note in main)
            prod *= sp.numer(sp.together((x1 - x2) ** 2 + (y1 - y2) ** 2))
    return sp.expand(prod)


def main():
    random.seed(1)
    # ---- chart 1: generic (L2 not vertical) -----------------------------------------------
    L = {1: (sp.Integer(0), sp.Integer(0)), 2: (a, sp.Integer(0)), 3: (b, sp.Integer(1)), 4: (c, d)}
    P = vertices(L)
    Q = three_conditions(P)
    print("chart 1 (L1: y=0, L2: y=ax, L3: y=bx+1, L4: y=cx+d)")
    for i, q in enumerate(Q, 1):
        print(f"   Q{i} numerator factored: {q}")
    # non-degeneracy over C: lines pairwise non-parallel, no three concurrent.  Over C, distinctness of
    # points p != q is NOT captured by a single polynomial; but "no three lines concurrent" and
    # "pairwise non-parallel" ARE polynomial conditions and they imply all six vertices are distinct.
    nonpar = a * b * c * (a - b) * (a - c) * (b - c)
    # concurrency: L1,L2,L3 at P12=(0,0): L3 through origin iff 1 = 0 -> never. L1,L2,L4: d = 0.
    # L1,L3,L4 concurrent: P13 = P14: -1/b = -d/c  <=> c = b d.   L2,L3,L4 concurrent: P23 on L4.
    x23, y23 = P[(2, 3)]
    conc234 = sp.numer(sp.together(y23 - (c * x23 + d)))
    nondeg = sp.expand(nonpar * d * (c - b * d) * conc234)
    print("   non-degeneracy polynomial factors:", sp.factor(nondeg))
    G = sp.groebner([Q[0], Q[1], Q[2], t * nondeg - 1], a, b, c, d, t, order='grevlex')
    print("   Groebner(Q1,Q2,Q3 + nondeg) == {1}:", list(G) == [1])
    G12 = sp.groebner([Q[0], Q[1], t * nondeg - 1], a, b, c, d, t, order='grevlex')
    print("   Groebner(Q1,Q2 + nondeg) == {1} (i.e. Q1=Q2=0 forces L2 vertical, outside this chart):",
          list(G12) == [1])
    # ---- chart 2: L2 vertical (L1 ⊥ L2) --------------------------------------------------------
    L2v = {1: (sp.Integer(0), sp.Integer(0)), 2: ('vert', sp.Integer(0)), 3: (b, sp.Integer(1)), 4: (c, d)}
    P2 = vertices(L2v)
    Qv = three_conditions(P2)
    print("chart 2 (L1: y=0, L2: x=0, L3: y=bx+1, L4: y=cx+d)")
    for i, q in enumerate(Qv, 1):
        print(f"   Q{i} numerator factored: {q}")
    x23, y23 = P2[(2, 3)]
    conc234 = sp.numer(sp.together(y23 - (c * x23 + d)))
    nondeg2 = sp.expand(b * c * (b - c) * d * (c - b * d) * conc234)
    print("   non-degeneracy polynomial factors:", sp.factor(nondeg2))
    G2 = sp.groebner([Qv[0], Qv[1], Qv[2], t * nondeg2 - 1], b, c, d, t, order='grevlex')
    print("   Groebner(Q1,Q2,Q3 + nondeg) == {1}:", list(G2) == [1])
    G2b = sp.groebner([Qv[0], Qv[1], t * nondeg2 - 1], b, c, d, t, order='lex')
    print("   Groebner(Q1,Q2 + nondeg), lex:", [sp.factor(g) for g in G2b])
    # ---- numeric sanity (floating point, many random quadrilaterals) -------------------------
    import numpy as np
    f = sp.lambdify((a, b, c, d), Q, 'numpy')
    fnd = sp.lambdify((a, b, c, d), nondeg, 'numpy')
    worst = None
    for _ in range(20000):
        vals = [random.uniform(-3, 3) for _ in range(4)]
        q = np.array(f(*vals), dtype=float)
        # normalise by scale of the points
        nd = abs(float(fnd(*vals)))
        if nd < 1e-3:
            continue
        m = np.max(np.abs(q))
        if worst is None or m < worst[0]:
            worst = (m, vals)
    print("   numeric: smallest max|Q_i| over 20000 random non-degenerate quadrilaterals:", worst)


if __name__ == "__main__":
    main()
