#!/usr/bin/env python3
"""analyze_C2.py -- exact realisability analysis of the last surviving n = 10 structure ("C-2"):
two disjoint 5-blocks A = {0..4}, B = {5..9} and twenty 4-blocks, each meeting A and B in a pair.

Moebius normal forms (A, B are two distinct circles on the sphere, no point of P on both):
 (i)   A, B disjoint  -> concentric circles |z| = 1 and |z| = r (r > 0, r != 1):  points e^{i a_j}, r e^{i b_k}.
       {a_i, a_j, b_k, b_l} concyclic  <=>  a_i + a_j == b_k + b_l (mod 2 pi)     [centre on the common bisector]
 (ii)  A, B tangent   -> parallel lines Im z = 0 and Im z = 1:  points x_j, t_k + i.
       concyclic  <=>  x_i + x_j = t_k + t_l                                       [centre on the common bisector]
 (iii) A, B meet in two points -> lines through 0 (real axis and a second line), points x_j, t_k e^{i theta}:
       concyclic  <=>  x_i x_j = t_k t_l   (signed coordinates; power of the origin), all coordinates nonzero.
So every realisation is a solution of a LINEAR system: over the torus (i), over R (ii), or (iii) over R for
log|.| together with signs.  We solve each exactly and list all solutions with the required distinctness.
Then, for every solution family, we compute the maximum number of blocks through a common point O of the
sphere (= the maximum number of lines after inverting in O), which gives the Euclidean circle count
42 - deg(O) of the corresponding planar realisations.
"""
import sys, json, itertools, math, cmath
from fractions import Fraction
import sympy as sp
import numpy as np

data = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'n10_table_C.json'))
rec = data['results'][int(sys.argv[2]) if len(sys.argv) > 2 else 2]
blocks = [tuple(b) for b in rec['blocks']]
A = [b for b in blocks if len(b) == 5][0]
B = [b for b in blocks if len(b) == 5][1]
fours = [b for b in blocks if len(b) == 4]
print("A =", A, " B =", B, " #4-blocks =", len(fours))
pairs = []
for b in fours:
    pa = tuple(sorted(set(b) & set(A))); pb = tuple(sorted(set(b) & set(B)))
    assert len(pa) == 2 and len(pb) == 2
    pairs.append((pa, pb))
# ---------------------------------------------------------------- integer matrix of the linear system
# unknown vector theta = (a_0..a_4, b_5..b_9); equation a_i + a_j - b_k - b_l = 0
idx = {p: i for i, p in enumerate(list(A) + list(B))}
M = sp.zeros(len(pairs), 10)
for r, (pa, pb) in enumerate(pairs):
    M[r, idx[pa[0]]] += 1; M[r, idx[pa[1]]] += 1; M[r, idx[pb[0]]] -= 1; M[r, idx[pb[1]]] -= 1
print("rank over Q of the 20 x 10 system:", M.rank())
# ---------------------------------------------------------------- case (ii): real solutions
ns = M.nullspace()
print(f"case (ii) tangent circles: real solution space of dimension {len(ns)}: basis {[list(v.T) for v in ns]}")
def distinct_ok(v):
    a = [v[i] for i in range(5)]; b = [v[i] for i in range(5, 10)]
    return len(set(a)) == 5 and len(set(b)) == 5
# any real solution is a combination of the basis; check whether some combination has distinct coordinates:
# the coordinates a_i - a_j are linear forms in the combination coefficients; if every basis vector has a_i = a_j
# for some fixed pair, no combination separates them.
if len(ns) <= 1:
    print("   -> only the constant solution (all points coincide): case (ii) impossible; case (iii) impossible too"
          " (|x_i| would all be equal, so at most 2 distinct signed coordinates on each line)")
else:
    lam = sp.symbols('l0:%d' % len(ns))
    v = sum((lam[i] * ns[i] for i in range(len(ns))), sp.zeros(10, 1))
    forms = [sp.simplify(v[i] - v[j]) for i, j in itertools.combinations(range(5), 2)] + \
            [sp.simplify(v[i] - v[j]) for i, j in itertools.combinations(range(5, 10), 2)]
    print("   difference forms:", forms)
    print("   some form identically zero:", any(f == 0 for f in forms))
# ---------------------------------------------------------------- case (i): torus solutions via Smith normal form
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy import ZZ
D = smith_normal_form(M, domain=ZZ)
divs = [D[i, i] for i in range(min(D.shape))]
print("Smith normal form diagonal:", divs)
# Solutions theta in (R/Z)^10 of M theta = 0 mod 1:  with M = U D V (U,V unimodular), phi = V theta,
# d_i phi_i in Z  =>  phi_i in (1/d_i) Z / Z  (finite) for d_i != 0,  phi_i free for d_i = 0.
# We compute V explicitly by column operations (sympy's smith_normal_form does not return transforms), so we
# instead enumerate solutions directly: the torsion part has exponent lcm(d_i); all torsion solutions are
# rational vectors with denominators dividing that exponent, and the continuous part is the real kernel.
exp_ = 1
for d in divs:
    if d != 0:
        exp_ = sp.ilcm(exp_, int(d))
print("exponent of the torsion part:", exp_, "; dimension of the continuous part:", 10 - sum(1 for d in divs if d != 0))
# enumerate torsion solutions modulo the continuous part: solve M theta = 0 mod 1 with theta in (1/e) Z^10,
# i.e. M (e theta) = 0 mod e over Z/eZ.  Since the continuous part is (presumably) the rotation a_i = b_k = c,
# fix a_0 = 0 to kill it.
e = int(exp_)
Mint = np.array(M.tolist(), dtype=np.int64)
sols = []
# solve over Z/eZ by brute force on the 9 free coordinates would be e^9; instead do a linear solve mod e via
# enumerating the kernel with a small search: e is small (expected 5 or 10).
import itertools as it
cont_dim = 10 - sum(1 for d in divs if d != 0)
assert cont_dim == 1, "continuous part larger than the rotation; handle separately"
for vec in it.product(range(e), repeat=9):
    th = np.array((0,) + vec, dtype=np.int64)
    if np.all((Mint @ th) % e == 0):
        sols.append(th)
print(f"torsion solutions with a_0 = 0 over Z/{e}Z: {len(sols)}")
good = []
for th in sols:
    a = [int(x) for x in th[:5]]; b = [int(x) for x in th[5:]]
    if len(set(a)) == 5 and len(set(b)) == 5:
        good.append((tuple(a), tuple(b)))
print(f"solutions with all five a's distinct and all five b's distinct (angles 2 pi * k / {e}): {len(good)}")
for g in good:
    print("   a =", g[0], " b =", g[1])
# ---------------------------------------------------------------- geometry of the concentric realisations
def blocks_of_config(pts, tol=1e-9):
    """all circles/lines through >= 3 of the points (numerically), as (frozenset of indices)."""
    n = len(pts)
    found = {}
    for i, j, k in itertools.combinations(range(n), 3):
        a, b, c = pts[i], pts[j], pts[k]
        d = 2 * ((a.real) * (b.imag - c.imag) + b.real * (c.imag - a.imag) + c.real * (a.imag - b.imag))
        if abs(d) < 1e-12:
            key = ('L', i, j, k)  # collinear: handled by cross-ratio below
        S = frozenset(q for q in range(n) if abs(cmath.phase(((pts[i] - pts[k]) * (pts[j] - pts[q])) / ((pts[i] - pts[q]) * (pts[j] - pts[k])) if q not in (i, j, k) else 1) % math.pi) < 1e-7 or (abs(cmath.phase(((pts[i] - pts[k]) * (pts[j] - pts[q])) / ((pts[i] - pts[q]) * (pts[j] - pts[k]))) % math.pi - math.pi) < 1e-7) or q in (i, j, k))
        found[S] = 1
    return list(found)

def circle_through(a, b, c):
    """(centre, radius) or None if collinear."""
    d = 2 * (a.real * (b.imag - c.imag) + b.real * (c.imag - a.imag) + c.real * (a.imag - b.imag))
    if abs(d) < 1e-13:
        return None
    a2, b2, c2 = abs(a) ** 2, abs(b) ** 2, abs(c) ** 2
    ux = (a2 * (b.imag - c.imag) + b2 * (c.imag - a.imag) + c2 * (a.imag - b.imag)) / d
    uy = (a2 * (c.real - b.real) + b2 * (a.real - c.real) + c2 * (b.real - a.real)) / d
    ctr = complex(ux, uy)
    return ctr, abs(a - ctr)

def max_degree(pts):
    """max number of blocks (circles/lines through >= 3 points) passing through a common point O not in P,
    found among pairwise intersections of blocks; returns (deg, O, ell_infinity = number of lines)."""
    n = len(pts)
    blocks = {}
    for i, j, k in itertools.combinations(range(n), 3):
        cc = circle_through(pts[i], pts[j], pts[k])
        if cc is None:
            key = ('L', round((pts[j] - pts[i]).imag / abs(pts[j] - pts[i]) * (1 if (pts[j] - pts[i]).real >= 0 else -1), 6),
                   round(((pts[i].conjugate() * (pts[j] - pts[i])).imag) / abs(pts[j] - pts[i]), 6))
        else:
            key = ('C', round(cc[0].real, 6), round(cc[0].imag, 6), round(cc[1], 6))
        blocks.setdefault(key, set()).update((i, j, k))
    circles = [k for k in blocks if k[0] == 'C']
    lines = [k for k in blocks if k[0] == 'L']
    # candidate points: intersections of pairs of circles
    cand = []
    for k1, k2 in itertools.combinations(circles, 2):
        c1 = complex(k1[1], k1[2]); r1 = k1[3]; c2 = complex(k2[1], k2[2]); r2 = k2[3]
        d = abs(c2 - c1)
        if d < 1e-9 or d > r1 + r2 + 1e-9 or d < abs(r1 - r2) - 1e-9:
            continue
        aa = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
        h2 = r1 * r1 - aa * aa
        h = math.sqrt(max(h2, 0.0))
        p = c1 + aa * (c2 - c1) / d
        for s in (1, -1):
            cand.append(p + s * h * 1j * (c2 - c1) / d)
    best = (len(lines), None)
    for O in cand:
        if min(abs(O - q) for q in pts) < 1e-6:
            continue
        deg = 0
        for k in circles:
            if abs(abs(O - complex(k[1], k[2])) - k[3]) < 1e-6:
                deg += 1
        if deg > best[0]:
            best = (deg, O)
    return best, len(circles), len(lines), len(blocks)

if good:
    a_ang, b_ang = good[0]
    print("\nnumerical scan of the concentric family for solution", good[0])
    results = []
    for r in np.linspace(0.05, 0.99, 95):
        pts = [cmath.exp(2j * math.pi * k / e) for k in a_ang] + [r * cmath.exp(2j * math.pi * k / e) for k in b_ang]
        (deg, O), nc, nl, nb = max_degree(pts)
        results.append((round(float(r), 3), deg, nb, nl))
    degs = {}
    for r, deg, nb, nl in results:
        degs.setdefault((deg, nb, nl), []).append(r)
    for k, v in sorted(degs.items()):
        print(f"   (max deg(O), #blocks, #lines) = {k}: r in {v[:6]}{'...' if len(v) > 6 else ''}")
