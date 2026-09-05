"""lemmas-audit: light but complete algebraic certificate of the Angle Lemma.

By the S4 symmetry of the statement (the three quadruples are the complements of the three
'opposite pairs' {P12,P34},{P13,P24},{P14,P23}, which S4 permutes), every configuration is, after
relabelling and a similarity, in one of two charts:
  chart A: L1: y=0 through P12=(0,0), and no line perpendicular to L1  ->  L_i: y = m_i x + c_i, c2=0.
  chart B: L1: y=0, L2: x=0 (a perpendicular pair, relabelled to be L1,L2), L3,L4: y = m_i x + c_i
           (they are neither horizontal nor vertical since they are not parallel to L1 or L2).
Chart A is the theory agent's chart; chart B is the one its certificate omits.
In each chart we factor the three concyclicity determinants, split off the factors that are
non-degeneracy conditions, and show the remaining 'essential' factors have no common zero
compatible with non-degeneracy (a Groebner basis of a 2-4 variable system).
"""
import itertools, sympy as sp

def concyc(pts):
    return sp.numer(sp.together(sp.Matrix([[x, y, x * x + y * y, 1] for (x, y) in pts]).det()))

def analyse_chart(name, lines, syms, nondeg):
    """lines: dict i -> (a,b,c) with a x + b y = c. nondeg: list of polynomials that must be nonzero."""
    def inter(i, j):
        a1, b1, c1 = lines[i]; a2, b2, c2 = lines[j]
        det = a1 * b2 - a2 * b1
        return sp.together((c1 * b2 - c2 * b1) / det), sp.together((a1 * c2 - a2 * c1) / det), det
    P = {}; dets = []
    for i, j in itertools.combinations(range(1, 5), 2):
        x, y, d = inter(i, j); P[(i, j)] = (x, y); dets.append(sp.factor(d))
    # non-degeneracy polynomials: dets (non-parallel) and vertex distinctness (numerators of dx, dy:
    # two vertices coincide iff both numerators vanish; since three lines concurrent <=> two vertices
    # coincide, we use the explicit 'concurrency' determinants instead, which are single polynomials).
    conc = []
    for i, j, k in itertools.combinations(range(1, 5), 3):
        M = sp.Matrix([list(lines[i]), list(lines[j]), list(lines[k])])
        conc.append(sp.factor(M.det()))
    nd_polys = set()
    for q in dets + conc + list(nondeg):
        for fac, _ in sp.factor_list(sp.expand(q))[1]:
            nd_polys.add(sp.Poly(fac, *syms).as_expr())
    Q = {1: concyc([P[(1, 2)], P[(3, 4)], P[(1, 3)], P[(2, 4)]]),
         2: concyc([P[(1, 2)], P[(3, 4)], P[(1, 4)], P[(2, 3)]]),
         3: concyc([P[(1, 3)], P[(2, 4)], P[(1, 4)], P[(2, 3)]])}
    ess = {}
    print(f"--- chart {name}: non-degeneracy factors: {sorted(map(str, nd_polys))}")
    for i, q in Q.items():
        fl = sp.factor_list(sp.expand(q))
        essential = []
        degenerate = []
        for fac, mult in fl[1]:
            fe = sp.Poly(fac, *syms).as_expr()
            if fe in nd_polys or sp.expand(-fe) in nd_polys:
                degenerate.append(fac)
            else:
                essential.append(fac)
        ess[i] = essential
        print(f"Q{i}: degenerate factors {degenerate}; essential factors {essential}")
    # common zeros of the essential factors: Groebner basis (lex), then inspect
    E = [sp.expand(e) for i in ess for e in ess[i]]
    G = sp.groebner(E, *syms, order="lex")
    print("Groebner basis (lex) of the essential factors:", list(G))
    # saturate by the non-degeneracy factors that actually matter: check each basis polynomial
    # against the nondeg set: if some product of non-degeneracy factors lies in the ideal we are done.
    t = sp.Symbol("t")
    prod = sp.Integer(1)
    for g in nd_polys:
        prod *= g
    G2 = sp.groebner(E + [t * prod - 1], *syms, t, order="grevlex")
    print("saturated by all non-degeneracy factors -> unit ideal:", list(G2) == [1])
    return list(G2) == [1]

if __name__ == "__main__":
    m2, m3, m4, c3, c4 = sp.symbols("m2 m3 m4 c3 c4")
    # chart A
    linesA = {1: (sp.Integer(0), sp.Integer(1), sp.Integer(0)), 2: (-m2, sp.Integer(1), sp.Integer(0)),
              3: (-m3, sp.Integer(1), c3), 4: (-m4, sp.Integer(1), c4)}
    okA = analyse_chart("A (theory's chart, no line perpendicular to L1)", linesA, [m2, m3, m4, c3, c4], [])
    # chart B
    linesB = {1: (sp.Integer(0), sp.Integer(1), sp.Integer(0)), 2: (sp.Integer(1), sp.Integer(0), sp.Integer(0)),
              3: (-m3, sp.Integer(1), c3), 4: (-m4, sp.Integer(1), c4)}
    okB = analyse_chart("B (L1 ⊥ L2, omitted by the theory's certificate)", linesB, [m3, m4, c3, c4], [])
    # In chart B, also show what Q1 = Q2 = 0 alone allows (should be consistent: a genuine family)
    print("ANGLE LEMMA PROVED IN BOTH CHARTS:", okA and okB)
