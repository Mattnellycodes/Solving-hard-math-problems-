"""Exact realisability analysis for structures whose blocks of size >= 5 are two DISJOINT 5-blocks S, R
(own derivation and code; the logic of n10-cpsat/v2/analyse_two5.py was re-derived, not copied).

Geometry (sphere/Möbius model).  S lies on a circle A, R on a circle B, A != B, and no point of P lies on both
(a point of P on A belongs to the block S).  A Möbius transformation brings (A, B) to a normal form:
 (i)   A, B disjoint  -> |z| = 1 and |z| = rho (rho > 0, rho != 1).  Points e^{2 pi i a}, rho e^{2 pi i b}.
       {a1, a2, b1, b2} concyclic/collinear  <=>  a1 + a2 == b1 + b2 (mod 1):  a circle through two points of
       |z| = 1 has its centre on the bisector of the chord (direction (a1+a2)/2 mod 1/2), likewise for the
       second chord; the two bisectors must coincide (the only circle centred at 0 through a point of A misses B),
       and conversely four points symmetric in a common axis form an isosceles trapezoid, hence are concyclic
       (or collinear).
 (ii)  A, B tangent      -> Im z = 0 and Im z = 1; points x, t + i; concyclic <=> x1 + x2 = t1 + t2.
 (iii) A, B secant       -> two lines through 0; points x (real) and t e^{i theta}; concyclic <=> x1 x2 = t1 t2
       (power of the point 0), so log|x1| + log|x2| = log|t1| + log|t2|.
Every 4-block is {s, s', r, r'} (a 4-block shares <= 2 points with S and with R), giving the row e_s + e_s' - e_r
- e_r' of an integer matrix M (rows = 4-blocks).  Realisations of the 4-block family are the solutions of
M theta = 0 in (R/Z)^10 (case i), R^10 (case ii; theta = positions), R^10 (case iii; theta = log|.|).
If ker_R M = span(1,...,1): cases (ii), (iii) force all five points of S to coincide (resp. to have equal modulus:
at most 2 distinct points) -> impossible; in case (i) every solution is a rotation of a torsion solution
theta ∈ (1/e) Z^10 (Smith normal form: M = U^-1 diag(d) V^-1, theta = V phi with d_i phi_i ∈ Z), e = d_max.
We enumerate the torsion solutions with distinct points on each circle, rotation fixed by theta_0 = 0.
Lines.  The Euclidean lines of the original picture are the blocks through the image O of infinity, a point of the
sphere not in P.  For each labelling and every line S ∈ L: O concyclic with three points p, q, r of S, i.e. the
cross-ratio (O, p; q, r) is real:  N conj(D) - conj(N) D = 0, N = (O - q)(p - r), D = (O - r)(p - q), with
u = O, v = conj(O) independent unknowns (complexification: no complex solution => no real one), and the case
O = infinity of the normal form (three points collinear).  Saturation by rho (rho^2 - 1) prod (u - z_p).
Groebner basis (1) in Q(zeta_e)[w, u, v, rho] for every labelling and both O-cases proves (F, L) unrealisable.
Usage: python3 analyse_two5_own.py post.json [indices...]   (structures with two disjoint 5-blocks)"""
import sys, json, itertools, math, time
import sympy as sp

# ---------------------------------------------------------------- Smith normal form with transforms (own code)
def smith(M):
    """U M V = D (diagonal, d_i | d_{i+1}), U, V unimodular; returns (D, U, V) as sympy integer matrices."""
    A = sp.Matrix(M); m, n = A.shape
    U = sp.eye(m); V = sp.eye(n)
    t = 0
    while t < min(m, n):
        # find a nonzero pivot with minimal absolute value in the submatrix
        best = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i, j] != 0 and (best is None or abs(A[i, j]) < abs(A[best[0], best[1]])):
                    best = (i, j)
        if best is None:
            break
        i, j = best
        A.row_swap(t, i); U.row_swap(t, i)
        A.col_swap(t, j); V.col_swap(t, j)
        while True:
            changed = False
            for i in range(t + 1, m):
                if A[i, t] != 0:
                    q = A[i, t] // A[t, t]
                    A[i, :] = A[i, :] - q * A[t, :]; U[i, :] = U[i, :] - q * U[t, :]
                    if A[i, t] != 0:
                        A.row_swap(t, i); U.row_swap(t, i); changed = True
            for j in range(t + 1, n):
                if A[t, j] != 0:
                    q = A[t, j] // A[t, t]
                    A[:, j] = A[:, j] - q * A[:, t]; V[:, j] = V[:, j] - q * V[:, t]
                    if A[t, j] != 0:
                        A.col_swap(t, j); V.col_swap(t, j); changed = True
            if not changed:
                # divisibility condition: make A[t,t] divide everything below-right
                bad = None
                for i in range(t + 1, m):
                    for j in range(t + 1, n):
                        if A[i, j] % A[t, t] != 0:
                            bad = (i, j); break
                    if bad: break
                if bad is None:
                    break
                A[t, :] = A[t, :] + A[bad[0], :]; U[t, :] = U[t, :] + U[bad[0], :]
        if A[t, t] < 0:
            A[t, :] = -A[t, :]; U[t, :] = -U[t, :]
        t += 1
    assert U * sp.Matrix(M) * V == A
    return A, U, V

def torsion_solutions(M):
    """all theta ∈ (Q/Z)^10 with M theta ∈ Z^m, as integer vectors v (theta = v/e), and e."""
    D, U, V = smith(M)
    m, n = D.shape
    d = [int(D[i, i]) for i in range(min(m, n)) if D[i, i] != 0]
    r = len(d)
    e = 1
    for x in d: e = math.lcm(e, x)
    sols = set()
    for ks in itertools.product(*[range(x) for x in d]):
        phi = [sp.Rational(k, x) for k, x in zip(ks, d)] + [0] * (n - r)
        theta = V * sp.Matrix(phi)
        v = tuple(int((theta[i] * e) % e) for i in range(n))
        sols.add(v)
    Mi = sp.Matrix(M)
    for v in sols:
        assert all(x % e == 0 for x in Mi * sp.Matrix(v))
    return sols, e, r

def analyse(rec, verbose=True):
    F = [frozenset(B) for B in rec['blocks']]; L = [frozenset(S) for S in rec['lines']]
    fives = [B for B in F if len(B) == 5]
    assert len(fives) == 2 and not (fives[0] & fives[1]) and all(len(B) in (4, 5) for B in F)
    S, R = sorted(fives[0]), sorted(fives[1])
    idx = {p: i for i, p in enumerate(S + R)}
    fours = [B for B in F if len(B) == 4]
    M = [[0] * 10 for _ in fours]
    for k, B in enumerate(fours):
        assert len(B & set(S)) == 2 and len(B & set(R)) == 2
        for p in B: M[k][idx[p]] = 1 if p in S else -1
    Ms = sp.Matrix(M)
    ker = Ms.nullspace()
    if verbose: print(f"   {len(fours)} four-blocks; rank {Ms.rank()}; real kernel dim {len(ker)}")
    if len(ker) != 1:
        print("   real kernel not 1-dimensional -> tangent/secant cases open -> NOT DECIDED"); return None
    assert list(ker[0].T) == [ker[0][0]] * 10
    sols, e, rank = torsion_solutions(M)
    # rotation normalisation: theta -> theta - theta_0 (rotation is a symmetry; M (1,..,1) = 0)
    norm = {tuple((x - s[0]) % e for x in s) for s in sols}
    good = sorted(s for s in norm if len(set(s[:5])) == 5 and len(set(s[5:])) == 5)
    if verbose:
        print(f"   torsion exponent e = {e}; torsion solutions {len(sols)}; with theta_0 = 0 and distinct points on each circle: {len(good)}")
        for s in good: print(f"      S angles {s[:5]}  R angles {s[5:]}   (units 2 pi / {e})")
    if not good:
        return True   # the 4-block family itself is unrealisable (any lines)
    zeta, u, v, rho, w = sp.symbols('zeta u v rho w')
    Phi = sp.cyclotomic_poly(e, zeta)
    def red(x): return sp.rem(sp.expand(x), Phi, zeta)
    all_dead = True
    for s in good:
        z = {}; zc = {}
        for p in S: z[p] = red(zeta ** s[idx[p]]); zc[p] = red(zeta ** ((e - s[idx[p]]) % e))
        for p in R: z[p] = red(rho * zeta ** s[idx[p]]); zc[p] = red(rho * zeta ** ((e - s[idx[p]]) % e))
        # sanity: 4-blocks concyclic
        for B in fours:
            a, b = sorted(B & set(S)); c, d = sorted(B & set(R))
            assert (s[idx[a]] + s[idx[b]] - s[idx[c]] - s[idx[d]]) % e == 0
        eqs = [Phi]
        for Sl in L:
            p, q, r_ = sorted(Sl)[:3]
            Nn = (u - z[q]) * (z[p] - z[r_]); Nc = (v - zc[q]) * (zc[p] - zc[r_])
            Dd = (u - z[r_]) * (z[p] - z[q]); Dc = (v - zc[r_]) * (zc[p] - zc[q])
            eqs.append(red(Nn * Dc - Nc * Dd))
        nd = rho * (rho - 1) * (rho + 1)
        for p in z: nd *= (u - z[p])
        eqs.append(red(sp.expand(nd * w - 1)))
        t0 = time.time()
        G = sp.groebner(eqs, w, u, v, rho, zeta, order='grevlex')
        fin = list(G.exprs) == [1]
        eqs2 = [Phi, red(sp.expand(rho * (rho - 1) * (rho + 1) * w - 1))]
        for Sl in L:
            p, q, r_ = sorted(Sl)[:3]
            d1 = z[q] - z[p]; d2 = z[r_] - z[p]; c1 = zc[q] - zc[p]; c2 = zc[r_] - zc[p]
            eqs2.append(red(d1 * c2 - c1 * d2))
        G2 = sp.groebner(eqs2, w, rho, zeta, order='grevlex')
        inf = list(G2.exprs) == [1]
        if verbose:
            print(f"      labelling S={s[:5]} R={s[5:]}: finite O impossible: {fin}; O = infinity impossible: {inf}  [{time.time()-t0:.1f}s]")
            if not fin: print("         finite-O basis:", list(G.exprs)[:6])
            if not inf: print("         infinity basis:", list(G2.exprs)[:6])
        all_dead &= (fin and inf)
    return all_dead

if __name__ == '__main__':
    data = json.load(open(sys.argv[1]))
    want = [int(x) for x in sys.argv[2:]]
    for rec in data:
        fives = [B for B in rec['blocks'] if len(B) == 5]
        if len(fives) != 2 or set(fives[0]) & set(fives[1]): continue
        if want and rec['index'] not in want: continue
        print(f"structure {rec['index']}: count {rec['count']}, l = {rec['l']}, lines {rec['lines']}")
        res = analyse(rec)
        print("   =>", "NOT REALISABLE (all labellings, finite and infinite O)" if res else "NOT DECIDED")
