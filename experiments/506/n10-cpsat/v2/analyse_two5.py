"""Exact realisability analysis for structures whose big blocks are two DISJOINT 5-blocks S, R (own code).

Geometry.  S lies on a circle A, R on a circle B (sphere model), A != B, no point of P on both.  A Möbius map
brings (A, B) to one of three normal forms:
  (i)  A, B disjoint:      |z| = 1 and |z| = rho (rho > 0, rho != 1); points e^{2 pi i a_j}, rho e^{2 pi i b_k};
       {a_i, a_j, b_k, b_l} concyclic-or-collinear  <=>  a_i + a_j == b_k + b_l (mod 1)
       (the circle/line through the four points is symmetric in the bisector through 0 of both chords);
  (ii) A, B tangent:       Im z = 0 and Im z = 1; points x_j, t_k + i;  concyclic  <=>  x_i + x_j = t_k + t_l;
  (iii) A, B secant:       two lines through 0; points x_j (real line), t_k e^{i theta};  concyclic  <=>
       x_i x_j = t_k t_l (power of 0), hence log|x_i| + log|x_j| = log|t_k| + log|t_l|.
So every 4-block {a_i,a_j,b_k,b_l} gives the same integer linear form a_i + a_j - b_k - b_l, and a realisation is a
solution of M theta = 0 over R/Z (i), over R (ii), or over R for the logs (iii).  If the real kernel of M is spanned by
the all-ones vector, (ii) and (iii) are impossible (points would coincide) and in (i) every solution is a rotation of
a torsion solution with denominators dividing e = exponent of the torsion group (Smith normal form); we enumerate
the torsion solutions mod e with distinct points on each circle.
Lines.  A line of the Euclidean picture is a block through the point O = image of infinity; O is a point of the
sphere not in P.  For each labelling and each line S in L we require O concyclic with three points of S:
cross-ratio(O, p; q, r) real, written with u = O and v = conj(O) as independent unknowns (complexification —
no complex solution implies no real one); plus the case O = infinity of the normal form (collinearity).  A
Groebner basis (1), after saturating rho (rho^2 - 1) and (O - points), proves (F, L) unrealisable.
Usage: python3 analyse_two5.py report.json   (output of postprocess.py; only surviving structures are analysed
unless --all is given)
"""
import sys, json, itertools, math, time, argparse
import numpy as np
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

ap = argparse.ArgumentParser(); ap.add_argument('file'); ap.add_argument('--all', action='store_true')
args = ap.parse_args()
data = json.load(open(args.file))

def kernel_mod(M, p):
    """all solutions of M v = 0 over Z/p (p prime), as tuples."""
    M = np.array(M, dtype=np.int64) % p; m, n = M.shape
    piv = []; row = 0
    for col in range(n):
        pr = next((r for r in range(row, m) if M[r, col] % p), None)
        if pr is None: continue
        M[[row, pr]] = M[[pr, row]]
        M[row] = (M[row] * pow(int(M[row, col]), -1, p)) % p
        for r in range(m):
            if r != row and M[r, col] % p:
                M[r] = (M[r] - M[r, col] * M[row]) % p
        piv.append(col); row += 1
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(n, dtype=np.int64); v[f] = 1
        for i, c in enumerate(piv): v[c] = (-M[i, f]) % p
        basis.append(v)
    sols = set()
    for coeffs in itertools.product(range(p), repeat=len(basis)):
        v = np.zeros(n, dtype=np.int64)
        for c, b in zip(coeffs, basis): v = (v + c * b) % p
        sols.add(tuple(int(x) for x in v))
    return sols

def torsion_solutions(M, e):
    """solutions of M theta == 0 (mod e) with theta in (Z/e)^10, via CRT over the prime powers of e (brute force
    lifting for prime powers)."""
    fac = sp.factorint(e)
    per = []
    for p, k in fac.items():
        q = p ** k
        sols = kernel_mod(M, p)
        # lift to p^k by brute force refinement
        for j in range(2, k + 1):
            qj = p ** j; new = set()
            for s in sols:
                for add in itertools.product(range(p), repeat=10):
                    v = tuple((s[i] + add[i] * p ** (j - 1)) % qj for i in range(10))
                    if all((sum(M[r][i] * v[i] for i in range(10))) % qj == 0 for r in range(len(M))):
                        new.add(v)
            sols = new
        per.append((q, sols))
    out = set()
    for combo in itertools.product(*[s for _, s in per]):
        v = []
        for i in range(10):
            val = 0
            for (q, _), s in zip(per, combo):
                # CRT
                Mq = e // q; inv = pow(Mq, -1, q)
                val += s[i] * Mq * inv
            v.append(val % e)
        out.add(tuple(v))
    assert all(all(sum(M[r][i] * v[i] for i in range(10)) % e == 0 for r in range(len(M))) for v in out)
    return out

u, v, rho, w, zeta = sp.symbols('u v rho w zeta')

def analyse(rec):
    F = [frozenset(B) for B in rec['blocks']]; L = [frozenset(S) for S in rec['lines']]
    fives = [B for B in F if len(B) == 5]
    assert len(fives) == 2 and not (fives[0] & fives[1]) and all(len(B) == 4 for B in F if len(B) != 5)
    S, R = sorted(fives[0]), sorted(fives[1])
    idx = {p: i for i, p in enumerate(S + R)}
    fours = [B for B in F if len(B) == 4]
    M = [[0] * 10 for _ in fours]
    for r_, B in enumerate(fours):
        assert len(B & set(S)) == 2 and len(B & set(R)) == 2
        for p in B: M[r_][idx[p]] += 1 if p in S else -1
    Ms = sp.Matrix(M)
    ker = Ms.nullspace()
    print(f"   {len(fours)} four-blocks; rank {Ms.rank()}; real kernel dimension {len(ker)}: {[list(k.T) for k in ker]}")
    if len(ker) != 1:
        print("   REAL KERNEL NOT 1-DIMENSIONAL: tangent/secant normal forms need separate treatment -> NOT DECIDED")
        return None
    D = smith_normal_form(Ms, domain=sp.ZZ)
    diag = [int(D[i, i]) for i in range(min(D.shape))]
    e = 1
    for d in diag:
        if d: e = math.lcm(e, abs(d))
    print(f"   Smith normal form diagonal {diag}; torsion exponent e = {e}")
    sols = torsion_solutions(M, e)
    good = sorted(s for s in sols if s[0] == 0 and len(set(s[:5])) == 5 and len(set(s[5:])) == 5)
    print(f"   torsion solutions mod {e}: {len(sols)}; with theta_0 = 0 and distinct points on each circle: {len(good)}")
    for s in good:
        print("      S angles (units 2pi/%d):" % e, s[:5], " R angles:", s[5:])
    # --- lines: O on all blocks of L, for each labelling
    Phi = sp.cyclotomic_poly(e, zeta)
    def red(expr): return sp.rem(sp.expand(expr), Phi, zeta)
    def conj(expr):
        # zeta -> zeta^(e-1), u <-> v, rho real
        return red(sp.expand(expr).subs({zeta: zeta ** (e - 1)}, simultaneous=True).subs({u: w}, simultaneous=True).subs({v: u}, simultaneous=True).subs({w: v}, simultaneous=True))
    all_impossible = True
    for s in good:
        pts = {}
        for p in S: pts[p] = red(zeta ** s[idx[p]])
        for p in R: pts[p] = red(rho * zeta ** s[idx[p]])
        # sanity: every 4-block concyclic in this labelling (angle condition)
        assert all((s[idx[a]] + s[idx[b]] - s[idx[c]] - s[idx[d]]) % e == 0 for B in fours for a, b in [sorted(B & set(S))] for c, d in [sorted(B & set(R))])
        eqs = [Phi]
        for Sl in L:
            p, q, r_ = [pts[x] for x in sorted(Sl)[:3]]
            Nn = (u - q) * (p - r_); Dd = (u - r_) * (p - q)
            eqs.append(red(Nn * conj(Dd) - conj(Nn) * Dd))
        nondeg = rho * (rho - 1) * (rho + 1)
        for p in pts: nondeg *= (u - pts[p])
        eqs.append(red(sp.expand(nondeg * w - 1)))
        t0 = time.time()
        G = sp.groebner(eqs, w, u, v, rho, zeta, order='grevlex')
        fin = list(G.exprs) == [1]
        # O = infinity of the normal form: lines are straight
        eqs2 = [Phi, red(sp.expand(rho * (rho - 1) * (rho + 1) * w - 1))]
        for Sl in L:
            p, q, r_ = [pts[x] for x in sorted(Sl)[:3]]
            d1 = q - p; d2 = r_ - p
            eqs2.append(red(d1 * conj(d2) - conj(d1) * d2))
        G2 = sp.groebner(eqs2, w, rho, zeta, order='grevlex')
        inf = list(G2.exprs) == [1]
        print(f"      labelling S={s[:5]} R={s[5:]}: finite O: basis (1) = {fin}; O = infinity: basis (1) = {inf}   [{time.time()-t0:.1f}s]")
        if not fin: print("         finite-O basis:", list(G.exprs)[:6])
        if not inf: print("         infinity basis:", list(G2.exprs)[:6])
        all_impossible &= (fin and inf)
    return all_impossible

for rec in data:
    fives = [B for B in rec['blocks'] if len(B) == 5]
    if len(fives) != 2 or set(fives[0]) & set(fives[1]): continue
    if not args.all and rec['report']['kills']: continue
    print(f"structure {rec['index']}: count {rec['count']}, lines {rec['lines']}")
    res = analyse(rec)
    print("   => NOT REALISABLE (all labellings, finite and infinite O)" if res else "   => NOT DECIDED")
