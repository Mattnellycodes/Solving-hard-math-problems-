"""Own exact realisability analysis for structures with two disjoint 5-blocks S = {0..4} (on a
circle A) and R = {5..9} (on a circle B != A; no point of P is on both), all 4-blocks {s,s',r,r'}.

Möbius normal forms of two distinct circles of the sphere (exhaustive):
 (i)  disjoint   -> |z| = 1 and |z| = rho (rho > 0, rho != 1).  {s,s',r,r'} concyclic  <=>
      alpha_s + alpha_s' = beta_r + beta_r' (mod 2 pi): the centre of the common circle lies on
      the perpendicular bisectors of both chords, which are lines through 0 (a line through s, s'
      is the limiting case); conversely reflection in that bisector maps r to r'.
 (ii) tangent    -> parallel lines y = 0, y = 1: concyclic <=> x_s + x_s' = t_r + t_r'.
 (iii) crossing  -> two lines through 0 (0 not in P): concyclic <=> x_s x_s' = t_r t_r'
      (power of the point 0), so log|x_s| + log|x_s'| = log|t_r| + log|t_r'|.
Hence in every case the 20 rows of the integer matrix M (+1 on s, s', -1 on r, r') annihilate the
parameter vector in an abelian group G: G = R/2piZ (i), R (ii), R x Z/2 (iii).  If the real kernel
of M is spanned by (1,...,1), (ii) and (iii) force the 5 points of S to coincide / have equal
modulus, impossible, and in (i) every solution is a rotation of a torsion solution, obtained from a
diagonalisation U M V = diag(d_i) (own implementation): angles/2pi = V y with y_i in (1/d_i) Z.

For the surviving structure(s) the lines are then tested exactly: the point O = image of the
original infinity must be co-circular with every line of L (O finite: unknowns z0 = O, w0 = conj O,
treated as independent, which only enlarges the solution set) or O = inf (lines are straight).
Everything is done in Q[zeta, rho, z0, w0, t] with the cyclotomic relation Phi_e(zeta) = 0 and the
saturation t rho (rho^2 - 1) = 1; a Gröbner basis {1} certifies that no complex solution exists.
"""
import json
import sys
import time
from fractions import Fraction
from itertools import combinations, product
from math import gcd
import sympy as sp
from lib import mask, bits, popcount, mobius_blocks, canon, check_structure

S = list(range(5))
R = list(range(5, 10))


# ---------------------------------------------------------------- integer diagonalisation
def diagonalise(M):
    """Unimodular U, V with U*M*V diagonal (integers).  Returns (U, D, V)."""
    A = [row[:] for row in M]
    m, n = len(A), len(A[0])
    U = [[int(i == j) for j in range(m)] for i in range(m)]
    V = [[int(i == j) for j in range(n)] for i in range(n)]

    def swap_rows(X, i, j):
        X[i], X[j] = X[j], X[i]

    def swap_cols(X, i, j):
        for row in X:
            row[i], row[j] = row[j], row[i]

    def add_row(X, i, j, c):  # row i += c * row j
        X[i] = [a + c * b for a, b in zip(X[i], X[j])]

    def add_col(X, i, j, c):  # col i += c * col j
        for row in X:
            row[i] += c * row[j]

    k = 0
    while k < min(m, n):
        # pivot: smallest nonzero |entry| in submatrix
        best = None
        for i in range(k, m):
            for j in range(k, n):
                if A[i][j] != 0 and (best is None or abs(A[i][j]) < abs(A[best[0]][best[1]])):
                    best = (i, j)
        if best is None:
            break
        i, j = best
        swap_rows(A, k, i); swap_rows(U, k, i)
        swap_cols(A, k, j); swap_cols(V, k, j)
        done = False
        while not done:
            done = True
            p = A[k][k]
            for i in range(k + 1, m):
                if A[i][k] != 0:
                    q = A[i][k] // p
                    add_row(A, i, k, -q); add_row(U, i, k, -q)
                    if A[i][k] != 0:
                        swap_rows(A, k, i); swap_rows(U, k, i)
                        done = False
                        break
            if not done:
                continue
            p = A[k][k]
            for j in range(k + 1, n):
                if A[k][j] != 0:
                    q = A[k][j] // p
                    add_col(A, j, k, -q); add_col(V, j, k, -q)
                    if A[k][j] != 0:
                        swap_cols(A, k, j); swap_cols(V, k, j)
                        done = False
                        break
        k += 1
    # verify
    Um, Mm, Vm = sp.Matrix(U), sp.Matrix(M), sp.Matrix(V)
    Dm = Um * Mm * Vm
    for i in range(m):
        for j in range(n):
            if i != j:
                assert Dm[i, j] == 0
    assert abs(Um.det()) == 1 and abs(Vm.det()) == 1
    return U, [Dm[i, i] for i in range(min(m, n))], V


def angle_labellings(F4):
    """All torsion solutions (angles/2pi as Fractions mod 1, rotation normalised alpha_0 = 0)
    with 5 distinct S-angles and 5 distinct R-angles.  Also returns rank and kernel info."""
    M = []
    for b in F4:
        row = [0] * 10
        for p in bits(b):
            row[p] = 1 if p < 5 else -1
        M.append(row)
    U, D, V = diagonalise(M)
    rank = sum(1 for d in D if d != 0)
    ker = sp.Matrix(M).nullspace()
    Vm = sp.Matrix(V)
    sols = []
    tors = [int(abs(d)) for d in D if d != 0]
    exponent = 1
    for d in tors:
        exponent = exponent * d // gcd(exponent, d)
    for ys in product(*[range(d) for d in tors]):
        y = [Fraction(yi, d) for yi, d in zip(ys, tors)] + [Fraction(0)] * (10 - rank)
        x = [sum(Fraction(int(Vm[i, j])) * y[j] for j in range(10)) for i in range(10)]
        x = [(xi - x[0]) % 1 for xi in x]
        if len({x[p] for p in S}) == 5 and len({x[p] for p in R}) == 5:
            sols.append(x)
    return {"rank": rank, "kernel_dim": len(ker), "kernel_is_constants": len(ker) == 1 and all(v == ker[0][0] for v in ker[0]),
            "torsion": tors, "exponent": exponent, "labellings": sols}


# ---------------------------------------------------------------- exact line analysis
def exact_kill(F, L, x, verbose=True):
    """Return (finite_ok, inf_ok): True means 'no solution' (Gröbner basis {1}) for that case."""
    dens = [xi.denominator for xi in x]
    e = 1
    for d in dens:
        e = e * d // gcd(e, d)
    zeta, rho, z0, w0, t = sp.symbols("zeta rho z0 w0 t")
    Phi = sp.Poly(sp.cyclotomic_poly(e, zeta), zeta)
    n_of = [int(xi * e) % e for xi in x]

    def red(expr):
        return sp.rem(sp.expand(expr), Phi.as_expr(), zeta)

    z = {}
    zb = {}
    for p in range(10):
        pw = zeta ** n_of[p]
        pwb = zeta ** ((e - n_of[p]) % e)
        if p < 5:
            z[p], zb[p] = pw, pwb
        else:
            z[p], zb[p] = rho * pw, rho * pwb
    eqs_fin, eqs_inf = [], []
    inf_possible = True
    for l in L:
        pts = bits(l)
        if l == mask(S):
            eqs_fin.append(z0 * w0 - 1)
            inf_possible = False
            continue
        if l == mask(R):
            eqs_fin.append(z0 * w0 - rho ** 2)
            inf_possible = False
            continue
        assert len(pts) in (3, 4)   # a 4-block used as a line: O on its circle (every triple)
        for (a, b, c) in combinations(pts, 3):
            rows = [[z[p] * zb[p], z[p], zb[p], 1] for p in (a, b, c)] + [[z0 * w0, z0, w0, 1]]
            eqs_fin.append(red(sp.Matrix(rows).det()))
            rows3 = [[z[p], zb[p], 1] for p in (a, b, c)]
            eqs_inf.append(red(sp.Matrix(rows3).det()))
    sat = t * rho * (rho ** 2 - 1) - 1
    G = sp.groebner(eqs_fin + [Phi.as_expr(), sat], z0, w0, t, rho, zeta, order="grevlex", domain="QQ")
    fin_ok = list(G.exprs) == [1]
    if inf_possible:
        Gi = sp.groebner(eqs_inf + [Phi.as_expr(), sat], t, rho, zeta, order="grevlex", domain="QQ")
        inf_ok = list(Gi.exprs) == [1]
    else:
        inf_ok = True
    if verbose:
        print(f"    labelling angles*{e} = {n_of}: finite O -> GB={'{1}' if fin_ok else G.exprs}; "
              f"O=inf -> {'impossible (5-line)' if not inf_possible else ('GB={1}' if inf_ok else Gi.exprs)}", flush=True)
    return fin_ok, inf_ok


def line_sets(F, need):
    covered = set()
    for b in F:
        for tr in combinations(bits(b), 3):
            covered.add(mask(tr))
    cands = list(F) + [mask(tr) for tr in combinations(range(10), 3) if mask(tr) not in covered]
    w = [popcount(c) - 1 if popcount(c) >= 4 else 3 for c in cands]
    pairs = [sp.binomial(popcount(c), 2) for c in cands]
    dF = [sum(sp.binomial(popcount(b) - 1, 2) for b in F if b >> p & 1) for p in range(10)]
    compat = [[popcount(a & b) <= 1 for b in cands] for a in cands]
    out = []
    m = len(cands)

    def rec(start, chosen, tot, q):
        if len(chosen) + (m - start) < need:
            return
        if len(chosen) >= need:
            out.append([cands[i] for i in chosen])
        for i in range(start, m):
            if any(not compat[i][j] for j in chosen) or tot + pairs[i] > 44:
                continue
            if any(q[p] + w[i] > 44 for p in bits(cands[i])):
                continue
            for p in bits(cands[i]):
                q[p] += w[i]
            chosen.append(i)
            rec(i + 1, chosen, tot + pairs[i], q)
            chosen.pop()
            for p in bits(cands[i]):
                q[p] -= w[i]

    rec(0, [], 0, list(dF))
    return out


if __name__ == "__main__":
    data = json.load(open(sys.argv[1]))  # stage-1 style file with "classes"
    mode = sys.argv[2] if len(sys.argv) > 2 else "survivor"
    t0 = time.time()
    n_real_labellings = 0
    summary = []
    for ci, cl in enumerate(data["classes"]):
        F = [mask(b) for b in cl["F"]]
        assert F[0] == mask(S) and F[1] == mask(R)
        F4 = F[2:]
        info = angle_labellings(F4)
        nl = len(info["labellings"])
        print(f"F-class {ci}: b4={len(F4)} rank={info['rank']} kernel_dim={info['kernel_dim']} constants={info['kernel_is_constants']} "
              f"torsion={info['torsion']} exponent={info['exponent']} distinct-point labellings={nl}", flush=True)
        assert info["kernel_is_constants"]
        summary.append((ci, nl))
        if nl == 0:
            continue
        n_real_labellings += 1
        need = 88 - (18 + 3 * len(F4))
        Ls = line_sets(F, need)
        print(f"  {len(Ls)} labelled line sets with l >= {need}", flush=True)
        for li, L in enumerate(Ls):
            check_structure(F, L)
            print(f"  line set {li}: {[bits(l) for l in L]}")
            for x in info["labellings"]:
                fin_ok, inf_ok = exact_kill(F, L, x)
                assert fin_ok and inf_ok, "NOT KILLED"
    print(f"F-classes with a distinct-point angle labelling: {n_real_labellings}; all their (labelled) line sets "
          f"and labellings killed exactly.  {time.time() - t0:.1f}s")
