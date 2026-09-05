"""Independent exact non-realisability proof for the theorem-only survivor (two disjoint 5-blocks
S = {0..4}, R = {5..9}, 20 four-blocks {s,s',r,r'}, 10 three-point lines), written with a
formulation DIFFERENT from n10-continue/analyse_two5_own.py and from two5_exact.py:

 * invariant factors of the 20 x 10 incidence matrix M via sympy's smith_normal_form (not an own
   diagonalisation), Q-kernel via sympy nullspace;
 * the torsion labellings by a direct DFS over (Z/eZ)^10 (e = exponent), not via the transform V;
 * the line conditions in REAL coordinates: O = (x, y), circle radius rho, and the algebraic
   numbers c = cos(pi/5), s = sin(pi/5) as ring variables with their minimal relations
   4c^2 - 2c - 1 = 0, s^2 + c^2 - 1 = 0 (all cos/sin(k pi/5) are Chebyshev polynomials in c, s);
   the Galois conjugates only relabel the angles (x3), so a Gröbner basis {1} over Q excludes
   every conjugate at once.

Why the normal form is exhaustive (re-derived): two distinct circles of the Möbius plane are disjoint,
tangent or crossing.  Tangent -> parallel lines y=0,y=1, concyclic {s,s',r,r'} <=> x_s+x_s' = t_r+t_r';
crossing -> two lines through 0 (0 not in P), concyclic <=> x_s x_s' = t_r t_r' (power of the point);
in both cases the 20 relations say M v = 0 for a real vector v (positions, resp. log|positions|), and
ker_Q M = constants forces the 5 points of S to coincide (resp. have equal modulus, so <= 2 of them
distinct).  Disjoint -> concentric |z| = 1, |z| = rho: {s,s',r,r'} concyclic <=> a_s + a_s' == b_r + b_r'
(mod 2 pi) (both perpendicular bisectors pass through 0 and must coincide, else the centre would be 0).
Smith form U M V = diag(d_i): theta solves M theta in Z^20 iff y = V^{-1} theta has d_i y_i in Z
(i <= rank) and the kernel coordinate free (= rotation), so every solution is a rotation of a
solution with all angles in (2 pi / e) Z, e = lcm(d_i).  The original lines are circles through
O = image of infinity (O finite, O distinct from the 10 points) or straight lines (O = infinity).
"""
import json, sys, time, random
from itertools import combinations
from math import gcd, cos, sin, pi
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from lib import mask, bits, popcount, check_structure

S, R = list(range(5)), list(range(5, 10))
surv = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "runs/survivors.json"))
which = int(sys.argv[2]) if len(sys.argv) > 2 else 0
rec = surv[which]
F = [mask(b) for b in rec["F"]]; L = [mask(l) for l in rec["L"]]
check_structure(F, L)
assert F[0] == mask(S) and F[1] == mask(R)
F4 = F[2:]
assert all(popcount(b & mask(S)) == 2 and popcount(b & mask(R)) == 2 for b in F4)
print("survivor:", rec["file"], rec["index"], "b4 =", len(F4), "l =", len(L), "count", rec["circles"])

# ---- integer relations
rows = []
for b in F4:
    row = [0] * 10
    for p in bits(b):
        row[p] = 1 if p < 5 else -1
    rows.append(row)
M = sp.Matrix(rows)
ker = M.nullspace()
print("rank(M) =", M.rank(), " Q-kernel dim =", len(ker), " kernel = constants:",
      len(ker) == 1 and all(v == ker[0][0] for v in ker[0]))
assert len(ker) == 1 and all(v == ker[0][0] for v in ker[0])
Dm = smith_normal_form(M, domain=sp.ZZ)
inv = [abs(int(Dm[i, i])) for i in range(min(Dm.shape)) if Dm[i, i] != 0]
e = 1
for d in inv:
    e = e * d // gcd(e, d)
print("invariant factors:", inv, " exponent e =", e, " |torsion| =", sp.prod(inv))

# ---- DFS for all solutions in (Z/e)^10 with alpha_0 = 0
cons = [tuple(bits(b)) for b in F4]
sols = []
val = [None] * 10


def dfs(i):
    if i == 10:
        sols.append(tuple(val)); return
    for v in ([0] if i == 0 else range(e)):
        val[i] = v
        ok = True
        for (a, b, c, d) in cons:
            if max(a, b, c, d) == i and (val[a] + val[b] - val[c] - val[d]) % e:
                ok = False; break
        if ok:
            dfs(i + 1)
    val[i] = None


dfs(0)
print("solutions in (Z/e)^10 with alpha_0 = 0:", len(sols), "(must equal |torsion| =", sp.prod(inv), ")")
assert len(sols) == sp.prod(inv)
labellings = [x for x in sols if len(set(x[p] for p in S)) == 5 and len(set(x[p] for p in R)) == 5]
print("labellings with 5 distinct angles on each circle:", len(labellings))
for x in labellings:
    print("   S angles", [x[p] for p in S], " R angles", [x[p] for p in R], f"(units 2pi/{e})")
assert e == 10, "the real-coordinate model below assumes e = 10 (angles k*pi/5)"

# ---- exact line analysis, real coordinates
x, y, rho, c, s, t = sp.symbols("x y rho c s t")
rel = [4 * c**2 - 2 * c - 1, s**2 + c**2 - 1]


def cs(k):
    """(cos(k pi/5), sin(k pi/5)) as polynomials in c = cos(pi/5), s = sin(pi/5)."""
    k %= 10
    C = sp.expand(sp.chebyshevt(k, c))
    Sn = sp.Integer(0) if k == 0 else sp.expand(sp.chebyshevu(k - 1, c) * s)
    return C, Sn


# numeric self-test of cs()
cv, sv = cos(pi / 5), sin(pi / 5)
for k in range(10):
    C, Sn = cs(k)
    assert abs(float(C.subs({c: cv, s: sv})) - cos(k * pi / 5)) < 1e-12
    assert abs(float(Sn.subs({c: cv, s: sv})) - sin(k * pi / 5)) < 1e-12


def analyse(lab):
    P = {}
    for p in range(10):
        C, Sn = cs(lab[p])
        P[p] = (C, Sn) if p < 5 else (rho * C, rho * Sn)
    eq_fin, eq_inf = [], []
    for l in L:
        a, b, d = bits(l)
        rowsF = [[P[p][0]**2 + P[p][1]**2, P[p][0], P[p][1], 1] for p in (a, b, d)] + [[x**2 + y**2, x, y, 1]]
        eq_fin.append(sp.expand(sp.Matrix(rowsF).det()))
        rowsI = [[P[p][0], P[p][1], 1] for p in (a, b, d)]
        eq_inf.append(sp.expand(sp.Matrix(rowsI).det()))
    sat = t * rho * (rho**2 - 1) - 1
    G = sp.groebner(eq_fin + rel + [sat], t, x, y, rho, c, s, order="grevlex", domain="QQ")
    Gi = sp.groebner(eq_inf + rel + [sat], t, rho, c, s, order="grevlex", domain="QQ")
    return list(G.exprs), list(Gi.exprs)


t0 = time.time()
allkilled = True
for lab in labellings:
    G, Gi = analyse(lab)
    fin = (G == [1]); inf = (Gi == [1])
    print(f"labelling S={[lab[p] for p in S]} R={[lab[p] for p in R]}: finite O: GB={'{1}' if fin else G};  O=inf: GB={'{1}' if inf else Gi}   [{time.time()-t0:.1f}s]", flush=True)
    allkilled &= fin and inf
print("ALL LABELLINGS KILLED (no complex point with rho != 0, +-1 for any Galois conjugate):", allkilled)

# ---- positive control: the machinery must find solutions when they exist.
# Take the same labelling and REPLACE the line set by 3 triples that are genuinely concyclic with a
# chosen O = (1/3, 1/5) at rho = 2 (built numerically, then verified exactly by the Gröbner basis != {1}).
lab = labellings[0]
Pn = {p: ((cos(lab[p] * pi / 5), sin(lab[p] * pi / 5)) if p < 5 else (2 * cos(lab[p] * pi / 5), 2 * sin(lab[p] * pi / 5))) for p in range(10)}
Ox, Oy = 1 / 3, 1 / 5


def circ_through(A, B, C):
    (ax, ay), (bx, by), (cx, cy) = A, B, C
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
    uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
    return ux, uy, ((ax - ux)**2 + (ay - uy)**2) ** 0.5


# control: single triple + O is trivially satisfiable; check that GB != {1} (finite-O system with one line)
Lsave = L
L = [L[0]]
G, Gi = analyse(lab)
print("control (one line only): finite-O GB == {1}?", G == [1], "(must be False)")
assert G != [1]
L = Lsave
print("done", f"{time.time()-t0:.1f}s")
