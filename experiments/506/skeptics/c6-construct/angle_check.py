"""Independent symbolic check of the Angle Lemma.
Lines L_i: y = m_i x + b_i (i=1..4). P_ij = L_i cap L_j.  Show:
 (1) the concyclicity determinant of each diagonal quadruple factors as
     (nondegeneracy factors) * (direction-only polynomial);
 (2) the three direction-only polynomials cannot vanish together unless two slopes coincide,
     via a Groebner basis of the ideal saturated by prod (m_i - m_j)  (Rabinowitsch trick).
"""
import sympy as sp
m1, m2, m3, m4, b1, b2, b3, b4, t = sp.symbols("m1 m2 m3 m4 b1 b2 b3 b4 t")
m = {1: m1, 2: m2, 3: m3, 4: m4}; b = {1: b1, 2: b2, 3: b3, 4: b4}


def P(i, j):
    x = (b[j] - b[i]) / (m[i] - m[j])
    return sp.together(x), sp.together(m[i] * x + b[i])


def cyc(pts):
    M = sp.Matrix([[x, y, x * x + y * y, 1] for x, y in pts])
    return sp.factor(sp.together(M.det()))


quads = {"{P12,P34,P13,P24}": [(1, 2), (3, 4), (1, 3), (2, 4)],
         "{P12,P34,P14,P23}": [(1, 2), (3, 4), (1, 4), (2, 3)],
         "{P13,P24,P14,P23}": [(1, 3), (2, 4), (1, 4), (2, 3)]}
conds = []
for name, q in quads.items():
    d = cyc([P(*ij) for ij in q])
    num, den = sp.fraction(d)
    print(name, "determinant numerator factors:", sp.factor(num))
    print("      denominator:", sp.factor(den))
    # pick the factor that involves no b's
    fac = [f for f, e in sp.factor_list(num)[1] if not (f.free_symbols & {b1, b2, b3, b4})]
    print("      b-free factors:", fac)
    nondeg = [f for f, e in sp.factor_list(num)[1] if (f.free_symbols & {b1, b2, b3, b4})]
    print("      b-dependent factors (degeneracies: coincident vertices):", nondeg)
    conds.append(fac)

# The direction-only conditions predicted by the angle argument: tan(t1+t4)=tan(t2+t3) etc.
pred = [(m1 + m4) * (1 - m2 * m3) - (m2 + m3) * (1 - m1 * m4),
        (m1 + m3) * (1 - m2 * m4) - (m2 + m4) * (1 - m1 * m3),
        (m1 + m2) * (1 - m3 * m4) - (m3 + m4) * (1 - m1 * m2)]
for c, p in zip(conds, pred):
    prodc = sp.Integer(1)
    for f in c:
        if f.free_symbols:
            prodc *= f
    q = sp.simplify(prodc / p)
    print("ratio (b-free factor)/(angle-predicted condition) =", q)

# Groebner: pred + t*prod(m_i - m_j) - 1  -> unit ideal?
disc = sp.Integer(1)
for i in range(1, 5):
    for j in range(i + 1, 5):
        disc *= (m[i] - m[j])
G = sp.groebner(pred + [t * disc - 1], t, m1, m2, m3, m4, order="grevlex")
print("Groebner basis of <A1,A2,A3, t*prod(mi-mj)-1>:", list(G))
# also the sums used in the angle argument
print("A1 - A2 =", sp.factor(pred[0] - pred[1]))
print("A1 - A3 =", sp.factor(pred[0] - pred[2]))
print("A1 + A2 =", sp.factor(pred[0] + pred[1]))
