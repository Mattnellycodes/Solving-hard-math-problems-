"""Exact analysis (auditor's own) of structures whose blocks of size >= 5 are two disjoint 5-blocks
S = {0..4}, R = {5..9}: F = {S, R} u {4-blocks {s,s',r,r'}} and a line set L.

Geometry.  In a realisation, S lies on a circle A and R on a circle B of the Möbius plane (a line
counts as a circle through oo), A != B, and no point of P is on both (else that circle would
contain >= 6 points).  A Möbius transformation T brings (A, B) to one of three normal forms:
  (i)   A, B disjoint  -> concentric circles |z| = 1, |z| = rho (rho > 0, rho != 1);
  (ii)  A, B tangent   -> parallel lines  y = 0, y = 1;
  (iii) A, B crossing  -> two lines through 0 (the crossing points go to 0 and oo).
In each case {s, s', r, r'} is concyclic (or collinear) iff
  (i)   alpha_s + alpha_s' == beta_r + beta_r'  (mod 2 pi)   [common perpendicular bisector of the
        chords ss', rr' through the centre; equal moduli excluded by rho != 1]
  (ii)  x_s + x_s' = t_r + t_r'                                 [isosceles trapezoid]
  (iii) x_s x_s' = t_r t_r'                                     [power of the point 0].
Let M be the b4 x 10 integer matrix with rows (+1 at s, s', -1 at r, r').  If ker M = span(1) then
(ii) forces all x_s equal and (iii) forces all |x_s| equal -- impossible for 5 distinct points -- so
only (i) remains, where the angle vector phi = theta / 2 pi satisfies M phi in Z^b4: modulo the
rotation (1,...,1) and modulo Z^10 the solutions are the finitely many torsion classes given by the
Smith normal form of M (snf.py).  A labelling is admissible only if the 5 angles of S are distinct
and the 5 angles of R are distinct.

Lines.  T(oo) = O.  Every line of L is a circle through oo, so its image is a generalised circle
through O: for each line, the 3 (or more) transformed points together with O lie on a common
circle-or-line; if O = oo the transformed points are collinear.  O is not an image of a point of P.
For each admissible labelling we test, exactly (Groebner bases over Q[.., z]/(Phi_N(z)) with z a
primitive N-th root of unity, N = lcm(e, 4), which also covers all Galois-conjugate labellings),
whether such an O exists for some rho with rho (rho^2 - 1) != 0:
  finite O = u + i v : det[(|p|^2, x_p, y_p, 1) for 3 points of the line; (u^2+v^2, u, v, 1)] = 0
  O = oo             : det[(x_p, y_p, 1) ...] = 0.
Groebner basis (1) after saturation by rho (rho^2 - 1) proves that no such O exists, hence that
(F, L) has no realisation.  If the basis is not (1), the solutions are printed for inspection
(spurious solutions have O equal to a point of P).

usage: python3 two5_exact.py fl_7.json [class indices...]      (default: all classes)
       python3 two5_exact.py --families f4_7.json              (F-classes only: labellings)
"""
import sys, json, itertools, math, time
from fractions import Fraction
import sympy as sp
from snf import snf, det as fdet
from common import N as NP, count_of

S = list(range(5)); R = list(range(5, 10))


def angle_matrix(F):
    rows = []
    for B in F:
        B = sorted(B)
        if len(B) != 4:
            continue
        s = [p for p in B if p < 5]; r = [p for p in B if p >= 5]
        assert len(s) == 2 and len(r) == 2, B
        row = [0] * 10
        row[s[0]] = row[s[1]] = 1; row[r[0]] = row[r[1]] = -1
        rows.append(row)
    return rows


def rank_frac(M):
    A = [[Fraction(x) for x in row] for row in M]
    rk = 0; m, n = len(A), len(A[0]); col = 0
    for col in range(n):
        p = next((r for r in range(rk, m) if A[r][col] != 0), None)
        if p is None:
            continue
        A[rk], A[p] = A[p], A[rk]
        for r in range(m):
            if r != rk and A[r][col] != 0:
                f = A[r][col] / A[rk][col]
                A[r] = [a - f * b for a, b in zip(A[r], A[rk])]
        rk += 1
    return rk


def labellings(F):
    """all admissible torsion labellings phi (tuples of Fractions mod 1, phi_0 = 0) and the exponent."""
    M = angle_matrix(F)
    rk = rank_frac(M)
    assert rk == 9, 'kernel larger than span(1): rank %d' % rk
    d, V = snf(M)
    d = [x for x in d if x != 0]
    assert len(d) == 9
    e = 1
    for x in d:
        e = e * x // math.gcd(e, x)
    labs = set()
    for ks in itertools.product(*[range(x) for x in d]):
        psi = [Fraction(k, x) for k, x in zip(ks, d)] + [Fraction(0)]
        phi = [sum(V[i][j] * psi[j] for j in range(10)) for i in range(10)]
        phi = [(p - phi[0]) % 1 for p in phi]          # rotate so that phi_0 = 0
        if len(set(phi[:5])) == 5 and len(set(phi[5:])) == 5:
            labs.add(tuple(phi))
    # sanity: every block satisfies the angle condition
    for phi in labs:
        for row in M:
            assert sum(a * b for a, b in zip(row, phi)) % 1 == 0
    return d, e, sorted(labs)


def numeric_check(F, phi, rho=1.7):
    """numerically verify that all 4-blocks are concyclic for this labelling (any rho)."""
    import cmath
    z = [cmath.exp(2j * math.pi * float(p)) * (1 if i < 5 else rho) for i, p in enumerate(phi)]
    for B in F:
        if len(B) == 4:
            a, b, c, d = [z[p] for p in B]
            cr = (a - c) * (b - d) / ((a - d) * (b - c))
            assert abs(cr.imag) < 1e-9, (B, cr)
    return True


def exact_test(F, L, phi, e, verbose=True):
    """Groebner test for O (finite and infinite).  Returns (basis_finite, basis_inf)."""
    Nn = e * 4 // math.gcd(e, 4)
    z, rho, u, v, t = sp.symbols('z rho u v t')
    Phi = sp.Poly(sp.cyclotomic_poly(Nn, z), z)
    ii = z ** (Nn // 4)                       # i = primitive 4th root

    def red(expr):
        return sp.rem(sp.Poly(sp.expand(expr), z, u, v, rho, t), sp.Poly(Phi.as_expr(), z, u, v, rho, t)).as_expr()

    coords = []
    for p in range(10):
        k = phi[p]
        assert (k * e).denominator == 1
        a = int(k * e) * (Nn // e) % Nn
        x = (z ** a + z ** ((Nn - a) % Nn)) / 2
        y = -ii * (z ** a - z ** ((Nn - a) % Nn)) / 2
        if p >= 5:
            x, y = rho * x, rho * y
        coords.append((sp.expand(x), sp.expand(y)))
    eqs_fin, eqs_inf = [], []
    for line in L:
        line = sorted(line)
        for tri in itertools.combinations(line, 3):
            rows = []
            for p in tri:
                x, y = coords[p]
                r2 = 1 if p < 5 else rho ** 2
                rows.append([r2, x, y, 1])
            Mf = sp.Matrix(rows + [[u ** 2 + v ** 2, u, v, 1]])
            eqs_fin.append(red(Mf.det()))
            Mi = sp.Matrix([[x, y, 1] for (x, y) in [coords[p] for p in tri]])
            eqs_inf.append(red(Mi.det()))
    sat = t * rho * (rho ** 2 - 1) - 1
    Gf = sp.groebner(eqs_fin + [Phi.as_expr(), sat], t, u, v, rho, z, order='lex', domain='QQ')
    Gi = sp.groebner(eqs_inf + [Phi.as_expr(), sat], t, rho, z, order='lex', domain='QQ')
    return Gf, Gi


def analyse_class(F, L, tag=''):
    d, e, labs = labellings(F)
    print('%s Smith diagonal %s exponent %d admissible labellings %d' % (tag, d, e, len(labs)), flush=True)
    results = []
    for phi in labs:
        numeric_check(F, phi)
        t0 = time.time()
        Gf, Gi = exact_test(F, L, phi, e)
        fin = list(Gf.exprs); inf = list(Gi.exprs)
        okf = (fin == [1]); oki = (inf == [1])
        print('   labelling %s: finite O -> %s ; O = oo -> %s  (%.0fs)' % (
            [str(p) for p in phi], '(1)' if okf else fin, '(1)' if oki else inf, time.time() - t0), flush=True)
        results.append((phi, okf, oki, fin, inf))
    return d, e, labs, results


if __name__ == '__main__':
    args = sys.argv[1:]
    if args[0] == '--families':
        data = json.load(open(args[1]))
        for i, cl in enumerate(data['classes']):
            F = cl['F']
            M = angle_matrix(F)
            rk = rank_frac(M)
            if rk != 9:
                print('F-class %d: rank %d (kernel larger than span(1)) -- needs separate treatment' % (i, rk))
                continue
            d, e, labs = labellings(F)
            print('F-class %d: b4=%d Smith diagonal %s exponent %d admissible labellings %d' % (i, cl['b4'], d, e, len(labs)))
        sys.exit()
    data = json.load(open(args[0]))
    idxs = [int(a) for a in args[1:]] or list(range(len(data['classes'])))
    for i in idxs:
        cl = data['classes'][i]
        F, L = cl['F'], cl['L']
        print('== class %d count=%d L=%s' % (i, count_of(F, L), L), flush=True)
        d, e, labs, res = analyse_class(F, L, tag='class %d:' % i)
        allk = all(okf and oki for _, okf, oki, _, _ in res)
        print('== class %d: %s' % (i, 'NO realisation (all labellings give (1) for finite and infinite O)' if allk else 'INSPECT'), flush=True)
