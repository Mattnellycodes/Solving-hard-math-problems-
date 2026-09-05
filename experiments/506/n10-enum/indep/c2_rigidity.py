"""Exact classification of the realisations of the surviving structure C-2 (candidate 4 of
n10_sg_t32.json): two disjoint 5-blocks A = {0..4}, B = {5..9} (circles with no common configuration
point) and 20 four-blocks {a,a',b,b'}.

Moebius normalisation of the two circles A, B (they share no point of P, so all points of P stay finite):
 (i)   A, B disjoint on the sphere  -> concentric: A = unit circle, B = circle |z| = rho, rho != 1.
       {a,a',b,b'} concyclic (circle or line)  <=>  alpha_a + alpha_a' == beta_b + beta_b' (mod 2 pi)
       (the chords aa', bb' of a circle meeting two concentric circles are both perpendicular to the
       line of centres; conversely the reflection in that line maps a<->a' and b to b').
 (ii)  A, B meet in two points  -> two lines through the origin: a = x e^{i0}, b = y e^{i gamma};
       concyclic <=> x x' = y y' (power of the origin), so log|x| + log|x'| = log|y| + log|y'|.
 (iii) A, B tangent -> two parallel lines y = 0, y = 1:  concyclic <=> x + x' = y + y'.
In all cases the incidences give the integer linear system M v = 0 (mod 2 pi in case (i), over R in
(ii), (iii)) with the same 20 x 10 matrix M.  This script computes: rank and real kernel of M (cases
(ii), (iii) need 5 distinct values on each line, impossible if the real kernel is the all-ones line),
the Smith normal form of M, and all solutions of M w == 0 (mod N) with distinct residues inside A and
inside B (N = exponent of the torsion part).  Output: c2_rigidity.json with the realisation types."""
import json, itertools, math, sys
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
import mo

data = json.load(open('n10_sg_t32.json'))['candidates']
C2 = [c for c in data if c['sizes'][:2] == [5, 5] and c['sizes'][2] == 4 and c['degrees'] == [9] * 10
      and not mo.miquel_violations(10, [mo.mask(b) for b in c['blocks']])]
assert len(C2) == 1, len(C2)
C2 = C2[0]
blocks = [tuple(b) for b in C2['blocks']]
A = set(blocks[0]); B = set(blocks[1]); assert len(A) == len(B) == 5 and not (A & B)
four = [b for b in blocks if len(b) == 4]
assert len(four) == 20 and all(len(set(b) & A) == 2 for b in four)
print("C-2 blocks:", blocks)
rows = []
for b in four:
    r = [0] * 10
    for p in b:
        r[p] = 1 if p in A else -1
    rows.append(r)
M = sp.Matrix(rows)
print("rank(M) over Q =", M.rank())
ker = M.nullspace()
print("real kernel basis:", [list(v.T) for v in ker])
assert M.rank() == 9 and len(ker) == 1 and all(x == ker[0][0] for x in ker[0])
print("=> cases (ii) and (iii) (intersecting / tangent A, B) impossible: all |x_a| resp. x_a equal, but A has 5 distinct points")
D = smith_normal_form(M, domain=sp.ZZ)
inv = [D[i, i] for i in range(min(D.shape)) if D[i, i] != 0]
print("Smith normal form invariant factors:", inv)
N = 1
for d in inv:
    N = sp.ilcm(N, abs(int(d)))
N = int(N)
print("torsion exponent N =", N)
# solve M w == 0 mod N by brute force over the kernel mod each prime factor (N squarefree assumed)
fac = sp.factorint(N)
assert all(e == 1 for e in fac.values()), fac
Mi = [[int(x) for x in row] for row in rows]

def kernel_mod_p(p):
    # Gaussian elimination over GF(p) on the 20x10 system; return all kernel vectors (list of tuples)
    A_ = [[x % p for x in row] for row in Mi]
    ncol = 10; piv = []; r = 0
    for c in range(ncol):
        pr = next((i for i in range(r, len(A_)) if A_[i][c] % p), None)
        if pr is None:
            continue
        A_[r], A_[pr] = A_[pr], A_[r]
        invv = pow(A_[r][c], -1, p)
        A_[r] = [(x * invv) % p for x in A_[r]]
        for i in range(len(A_)):
            if i != r and A_[i][c] % p:
                f = A_[i][c]
                A_[i] = [(x - f * y) % p for x, y in zip(A_[i], A_[r])]
        piv.append(c); r += 1
    free = [c for c in range(ncol) if c not in piv]
    sols = []
    for vals in itertools.product(range(p), repeat=len(free)):
        w = [0] * ncol
        for c, v in zip(free, vals):
            w[c] = v
        for i, c in enumerate(piv):
            w[c] = (-sum(A_[i][f] * w[f] for f in free)) % p
        sols.append(tuple(w))
    return sols

kers = {p: kernel_mod_p(p) for p in fac}
print("kernel sizes:", {p: len(k) for p, k in kers.items()})
# CRT combine
sols = []
for combo in itertools.product(*[kers[p] for p in fac]):
    w = []
    for i in range(10):
        x = 0
        for p, vec in zip(fac, combo):
            # CRT: find x mod N with x == vec[i] mod p
            Np = N // p
            x += vec[i] * Np * pow(Np, -1, p)
        w.append(x % N)
    assert all(sum(m * wi for m, wi in zip(row, w)) % N == 0 for row in Mi)
    sols.append(tuple(w))
print("solutions of M w == 0 (mod N):", len(sols))
Al = sorted(A); Bl = sorted(B)
good = [w for w in sols if len({w[a] for a in Al}) == 5 and len({w[b] for b in Bl}) == 5]
print("with distinct residues within A and within B:", len(good))
types = []
for w in good:
    if w[Al[0]] != 0:
        continue   # use the rotation to fix a_0 = 0
    ang = {p: 360.0 * w[p] / N for p in range(10)}
    regA = sorted(w[a] * (360 / N) % 360 for a in Al); regB = sorted(w[b] * (360 / N) % 360 for b in Bl)
    isregA = all(abs((regA[i + 1] - regA[i]) - 72) < 1e-9 for i in range(4))
    isregB = all(abs((regB[i + 1] - regB[i]) - 72) < 1e-9 for i in range(4))
    delta = (regB[0] - regA[0]) % 72
    types.append({'w': list(w), 'anglesA': [ang[a] for a in Al], 'anglesB': [ang[b] for b in Bl],
                  'A_regular_pentagon': isregA, 'B_regular_pentagon': isregB, 'B_offset_deg': delta})
    print(f"  solution (a_0 = 0): A angles {[ang[a] for a in Al]}  B angles {[ang[b] for b in Bl]} | A regular: {isregA} | B regular: {isregB} | B offset mod 72: {delta}")
json.dump({'blocks': [list(b) for b in blocks], 'A': Al, 'B': Bl, 'N': N, 'invariant_factors': [int(d) for d in inv],
           'types': types}, open('c2_rigidity.json', 'w'), indent=1)
print("CONCLUSION: every realisation of C-2 (points distinct, the 22 blocks concyclic) is Moebius-equivalent to one of the",
      len(types), "concentric configurations above, with an arbitrary radius ratio rho != 1 (and an overall rotation).")
