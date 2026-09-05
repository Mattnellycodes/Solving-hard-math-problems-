#!/usr/bin/env python3
"""Route C: the auditor's INDEPENDENT rigorous enclosure of

    G(k) := (-1)^k F(p_k#) = sum_{d | p_k#} (-1)^{omega(d)} H_d ,   k = 2..200,

written from scratch, with a different method mix than the certificate's routes (a)/(b)/checker:
  * arithmetic: exact python Fractions for every algebraic quantity; for the two transcendental inputs
    (ln p and gamma) my own interval class on top of GNU MPFR (gmpy2) with directed rounding (RoundDown /
    RoundUp contexts); every exact->mpfr conversion is checked against the exact rational.
  * ln p by the telescoping chain ln p = sum_{i=2}^{p} 2 atanh(1/(2i-1)), atanh(1/n) = sum_j 1/((2j+1) n^{2j+1})
    (only reciprocals of exact integers are rounded); cross-checked against MPFR's correctly rounded log.
  * gamma by TWO methods: (G1) gamma = Ein(x) - ln x - E_1(x), x = 256, Ein by its exact alternating power
    series, 0 < E_1(x) < e^{-x}/x < 2.7^{-x}/x  (A&S 5.1.11; no Euler-Maclaurin involved);
    (G2) gamma = H_N - ln N - A_m(N) - R_m(N), N = 2^16, m = 10 (lemma L1), H_N in interval arithmetic.
  * parameters (m, D) = (10, 3000) instead of (4, 10^5) / (6, 10^4); tail constant |B_22|/22.
  * divisor sets enumerated by factoring every integer d <= D (not by DFS over primes), with a brute-force
    recount for several k.
Mathematics (my own derivation, see README.md): for squarefree n = p_k#, k >= 2,
  G(k) = sum_s c_s P_s(k) + sum_{d<=D, d|n} (-1)^omega (H_d - ln d - gamma - A_m(d)) + E,
  |E| <= C (prod_{p<=p_k}(1+p^{-22}) - sum_{d<=D, d|n} d^{-22}),  C = |B_22|/22,
  c_1 = 1/2, c_{2j} = -B_{2j}/(2j), P_s(k) = prod_{i<=k}(1 - p_i^{-s}), A_m(d) = sum_s c_s d^{-s}.
"""
import sys, json, time, warnings, math
from fractions import Fraction
import gmpy2
from gmpy2 import mpfr, mpq, mpz
import sympy
sys.set_int_max_str_digits(0)
warnings.simplefilter("ignore", DeprecationWarning)
HERE = "/home/user/Solving-hard-math-problems-/experiments/A067857/audit-certificate"
CERT = "/home/user/Solving-hard-math-problems-/experiments/A067857/certificate"
M, D, KMAX, PREC = 10, 3000, 200, 420
S_TAIL = 2 * M + 2
t0 = time.time()

# ------------------------------------------------------------------ MPFR directed-rounding intervals
CD = gmpy2.context(precision=PREC, round=gmpy2.RoundDown)
CU = gmpy2.context(precision=PREC, round=gmpy2.RoundUp)
def frac(x):
    return Fraction(*x.as_integer_ratio())
def down_q(q):
    with gmpy2.local_context(CD):
        v = mpfr(mpq(q.numerator, q.denominator))
    assert frac(v) <= q
    return v
def up_q(q):
    with gmpy2.local_context(CU):
        v = mpfr(mpq(q.numerator, q.denominator))
    assert frac(v) >= q
    return v
def exact_small(c):
    v = mpfr(c)                       # default 53-bit context: exact for |c| < 2^53
    assert int(v) == c and abs(c) < 2 ** 53
    return v
class I:
    __slots__ = ("lo", "hi")
    def __init__(self, lo, hi):
        assert lo <= hi
        self.lo, self.hi = lo, hi
    @staticmethod
    def of_frac(q):
        return I(down_q(q), up_q(q))
    @staticmethod
    def of_frac_pair(lo, hi):
        return I(down_q(lo), up_q(hi))
    def __add__(a, b):
        return I(CD.add(a.lo, b.lo), CU.add(a.hi, b.hi))
    def __sub__(a, b):
        return I(CD.sub(a.lo, b.hi), CU.sub(a.hi, b.lo))
    def scale(a, c):                  # multiply by a small python int c (no unary minus: gmpy2 would apply the
        cm = exact_small(c)           # 53-bit GLOBAL context to "-x"; found by the containment check at k = 7)
        if c >= 0:
            return I(CD.mul(a.lo, cm), CU.mul(a.hi, cm))
        return I(CD.mul(a.hi, cm), CU.mul(a.lo, cm))
    def to_frac(a):
        return frac(a.lo), frac(a.hi)
    def width(a):
        return float(CU.sub(a.hi, a.lo))
    def contains(a, q):
        lo, hi = a.to_frac(); return lo <= q <= hi
def intersects(p, q):
    return p[0] <= q[1] and q[0] <= p[1]
ZERO = I(mpfr(0), mpfr(0))

# ------------------------------------------------------------------ primes (own sieve) + cross-check with sympy
LIM = 1300
sv = bytearray([1]) * (LIM + 1); sv[0] = sv[1] = 0
for i in range(2, int(LIM ** 0.5) + 1):
    if sv[i]:
        sv[i * i::i] = bytearray(len(sv[i * i::i]))
primes = [i for i in range(LIM + 1) if sv[i]][:KMAX]
assert len(primes) == KMAX and primes == list(sympy.primerange(2, 1224)) and sympy.prime(200) == 1223
assert primes[90] == 467 and primes[91] == 479 and primes[92] == 487
pindex = {p: i for i, p in enumerate(primes)}          # 0-based index

# ------------------------------------------------------------------ Bernoulli numbers by Akiyama-Tanigawa; cross-check sympy
def bernoulli_AT(nmax):
    B = []
    for n in range(nmax + 1):
        A = [Fraction(0)] * (n + 1)
        for mm in range(n + 1):
            A[mm] = Fraction(1, mm + 1)
            for j in range(mm, 0, -1):
                A[j - 1] = j * (A[j - 1] - A[j])
        B.append(A[0])            # Akiyama-Tanigawa gives B_n with B_1 = +1/2
    return B
Bl = bernoulli_AT(S_TAIL)
B = {2 * j: Bl[2 * j] for j in range(1, M + 2)}
for s, v in B.items():
    sb = sympy.bernoulli(s); assert Fraction(int(sb.p), int(sb.q)) == v
assert B[2] == Fraction(1, 6) and B[22] == Fraction(854513, 138) and B[20] == Fraction(-174611, 330)
coef = {1: Fraction(1, 2)}
for j in range(1, M + 1):
    coef[2 * j] = -B[2 * j] / (2 * j)
S_LIST = sorted(coef)
C_TAIL = abs(B[S_TAIL]) / S_TAIL
def A_m(d):
    return sum(c * Fraction(1, d ** s) for s, c in coef.items())
print("Bernoulli/coefficients ok; C_tail = |B_22|/22 = %s = %.4f   [%.1fs]" % (C_TAIL, float(C_TAIL), time.time() - t0))

# ------------------------------------------------------------------ ln n for n <= p_200 by the atanh chain
def atanh_recip(n):
    """Enclosure of atanh(1/n) = sum_{j>=0} 1/((2j+1) n^{2j+1}), integer n >= 3."""
    s = ZERO
    j, npow = 0, n
    limit = 1 << (PREC + 40)
    while True:
        s = s + I.of_frac(Fraction(1, (2 * j + 1) * npow))
        j += 1; npow *= n * n
        if (2 * j + 1) * npow > limit:
            break
    tail_hi = Fraction(n * n, n * n - 1) / ((2 * j + 1) * npow)   # sum_{i>=j} <= first term * 1/(1 - n^-2)
    return s + I(mpfr(0), up_q(tail_hi))
LNINT = [None, ZERO]
for i in range(2, primes[-1] + 1):
    LNINT.append(LNINT[i - 1] + atanh_recip(2 * i - 1).scale(2))
LN = {p: LNINT[p] for p in primes}
LN2 = LN[2]
# cross-check 1: MPFR correctly-rounded log (independent algorithm) must intersect
for p in primes:
    mine = LN[p].to_frac()
    lib = (frac(CD.log(exact_small(p))), frac(CU.log(exact_small(p))))
    assert intersects(mine, lib), p
print("ln p (chain of %d atanh series) widths <= %.2e, all consistent with MPFR log   [%.1fs]"
      % (primes[-1] - 1, max(LN[p].width() for p in primes), time.time() - t0))

# ------------------------------------------------------------------ gamma, method G1 (Ein / E_1), no Euler-Maclaurin
X = 256
Sk, term, k = Fraction(0), Fraction(1), 0
Kmax = 1200
while True:
    k += 1
    term = term * X / k                       # X^k / k!
    Sk += (-1) ** (k - 1) * term / k          # X^k / (k k!)
    if k > X and term / (k + 1) * X / (k + 1) < Fraction(1, 10 ** 130):
        break
tnext = term * X / (k + 1) / (k + 1)          # first omitted term (alternating, decreasing for k > X): |remainder| <= tnext
EB = Fraction(10, 27) ** X / X                # E_1(256) < e^{-256}/256 < 2.7^{-256}/256  (e > 2.7)
GAMMA_G1 = I.of_frac_pair(Sk - tnext - 8 * frac(LN2.hi) - EB, Sk + tnext - 8 * frac(LN2.lo))
# ------------------------------------------------------------------ gamma, method G2 (lemma L1 at N = 2^16, m = 10)
N = 1 << 16
HN = ZERO
for i in range(1, N + 1):
    HN = HN + I.of_frac(Fraction(1, i))
RN = C_TAIL / Fraction(N) ** S_TAIL
GAMMA_G2 = HN - LN2.scale(16) - I.of_frac(A_m(N)) - I.of_frac_pair(-RN, RN)
g1, g2 = GAMMA_G1.to_frac(), GAMMA_G2.to_frac()
assert intersects(g1, g2), "gamma enclosures G1 (Ein) and G2 (Euler-Maclaurin) disagree"
GAMMA = I.of_frac_pair(max(g1[0], g2[0]), min(g1[1], g2[1]))        # intersection (both are rigorous)
KNOWN = Fraction("0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329174674951463144724980708248096050401448654283622417399764492353625350033374293733773767394279259525824709491600873520394816567085323315177661152862119950150798479374508570574002992135478614669402960432542151905877553526733139925401296742051375413954911168510280798423487758720503843109399736137255306088933126760017247953783675927135157722610273492913940798430103417771778088154957066107501016191663340152278935867965497252036212879226555953669628176388623")
assert GAMMA.contains(KNOWN)
print("gamma: G1 (Ein, no E-M) width %.2e, G2 (E-M at 2^16, m=10) width %.2e, consistent; intersection width %.2e; contains published digits   [%.1fs]"
      % (GAMMA_G1.width(), GAMMA_G2.width(), GAMMA.width(), time.time() - t0))

# cross-check 2: certificate constants must intersect mine
for route in ("a", "b"):
    J = json.load(open("%s/cert_route_%s.json" % (CERT, route)))
    cg = (Fraction(J["constants"]["gamma"][0]), Fraction(J["constants"]["gamma"][1]))
    assert intersects(cg, GAMMA.to_frac()), route
    for p in primes:
        cl = (Fraction(J["constants"]["ln_p"][str(p)][0]), Fraction(J["constants"]["ln_p"][str(p)][1]))
        assert intersects(cl, LN[p].to_frac()), (route, p)
print("certificate constants (gamma, ln p for 200 primes, routes a and b) all intersect my enclosures")

# ------------------------------------------------------------------ exact Euler products, incremental in k
P = {s: [Fraction(1)] for s in S_LIST}         # P[s][k] = prod_{i<=k} (1 - p_i^{-s})
PT = [Fraction(1)]                               # prod_{i<=k} (1 + p_i^{-22})
for k in range(1, KMAX + 1):
    p = primes[k - 1]
    for s in S_LIST:
        P[s].append(P[s][-1] * (1 - Fraction(1, p ** s)))
    PT.append(PT[-1] * (1 + Fraction(1, p ** S_TAIL)))
print("exact Euler products done   [%.1fs]" % (time.time() - t0))

# ------------------------------------------------------------------ exact H_d, d <= D
H = [Fraction(0)]
for i in range(1, D + 1):
    H.append(H[-1] + Fraction(1, i))
hs = sympy.harmonic(D); assert Fraction(int(hs.p), int(hs.q)) == H[D]

# ------------------------------------------------------------------ divisor data by factoring each d <= D
def factor(n):
    f, q = [], 2
    while q * q <= n:
        while n % q == 0:
            f.append(q); n //= q
        q += 1
    if n > 1:
        f.append(n)
    return f
info = {}                                          # d -> (omega, gpf_index(1-based), prime indices)
for d in range(1, D + 1):
    f = factor(d)
    if len(set(f)) != len(f):
        continue                                    # not squarefree
    if f and f[-1] > primes[-1]:
        continue                                    # prime factor beyond p_200: does not divide p_200#
    idx = [pindex[q] for q in f]
    info[d] = (len(f), (max(idx) + 1) if idx else 0, idx)
assert info[1] == (0, 0, [])
V = {d: (-1) ** w * (H[d] - A_m(d)) for d, (w, g, idx) in info.items()}
W = {d: Fraction(1, d ** S_TAIL) for d in info}
print("%d squarefree p_200-smooth d <= %d   [%.1fs]" % (len(info), D, time.time() - t0))

# ------------------------------------------------------------------ per-k accumulation (prefix over gpf index)
by_g = {}
for d, (w, g, idx) in info.items():
    by_g.setdefault(g, []).append(d)
vsum, wsum, c0 = Fraction(0), Fraction(0), 0
cp = [0] * KMAX
ndiv = 0
def absorb(g):
    global vsum, wsum, c0, ndiv
    for d in by_g.get(g, []):
        w, _, idx = info[d]
        sg = (-1) ** w
        vsum += V[d]; wsum += W[d]; c0 += sg; ndiv += 1
        for i in idx:
            cp[i] += sg
absorb(0)
certA = {c["k"]: c for c in json.load(open("%s/cert_route_a.json" % CERT))["certificates"]}
certB = {c["k"]: c for c in json.load(open("%s/cert_route_b.json" % CERT))["certificates"]}
def dec(q, digits, up):
    n = -((-q.numerator * 10 ** digits) // q.denominator) if up else (q.numerator * 10 ** digits) // q.denominator
    s = str(abs(n)).rjust(digits + 1, "0")
    return ("-" if n < 0 else "") + s[:-digits] + "." + s[-digits:]
results = []
first_neg = None
for k in range(1, KMAX + 1):
    absorb(k)
    if k < 2:
        continue
    # brute-force recount of the divisor set for a few k (independent of the prefix bookkeeping)
    if k in (2, 3, 5, 10, 91, 92, 93, 200):
        pk = primes[k - 1]
        Sk = [d for d in range(1, D + 1) if d in info and (d == 1 or max(factor(d)) <= pk)]
        assert len(Sk) == ndiv
        assert sum(V[d] for d in Sk) == vsum and sum(W[d] for d in Sk) == wsum
        assert sum((-1) ** info[d][0] for d in Sk) == c0
        for i in range(k):
            assert sum((-1) ** info[d][0] for d in Sk if i in info[d][2]) == cp[i]
        assert all(cp[i] == 0 for i in range(k, KMAX))
    closed = sum(coef[s] * P[s][k] for s in S_LIST)
    Q = closed + vsum
    tail = PT[k] - wsum
    assert tail >= 0
    T = C_TAIL * tail
    Xi = GAMMA.scale(c0)
    for i in range(k):
        if cp[i]:
            Xi = Xi + LN[primes[i]].scale(cp[i])
    Xlo, Xhi = Xi.to_frac()
    Glo, Ghi = Q - Xhi - T, Q - Xlo + T
    sgn = -1 if Ghi < 0 else (1 if Glo > 0 else 0)
    if sgn < 0 and first_neg is None:
        first_neg = k
    a_sign = (-1) ** k * sgn
    rec = {"k": k, "p_k": primes[k - 1], "G_lo": dec(Glo, 100, False), "G_hi": dec(Ghi, 100, True),
           "width": float(Ghi - Glo), "tail_T": float(T), "sign_G": sgn, "sign_a": a_sign,
           "num_divisors_exact": ndiv, "c0": c0}
    for name, cert in (("a", certA), ("b", certB)):
        c = cert[k]
        clo, chi = Fraction(c["G_enclosure"][0]), Fraction(c["G_enclosure"][1])
        rec["inside_route_" + name] = (clo <= Glo and Ghi <= chi)
        rec["intersects_route_" + name] = intersects((Glo, Ghi), (clo, chi))
        rec["sign_agrees_route_" + name] = (c["sign_G"] == sgn and c["sign_a"] == a_sign)
    results.append(rec)
    if k in (2, 3, 4, 5, 6, 7, 50, 88, 89, 90, 91, 92, 93, 94, 95, 96, 100, 150, 200):
        print("k=%3d p_k=%4d  G in [%s, %s]  width %.1e  T=%.1e  sign %+d  inside a/b: %s/%s"
              % (k, primes[k - 1], dec(Glo, 22, False), dec(Ghi, 22, True), float(Ghi - Glo), float(T), sgn,
                 rec["inside_route_a"], rec["inside_route_b"]))
print("first k with G(k) < 0:", first_neg, "   [%.1fs]" % (time.time() - t0))
from mpmath import mp, harmonic as mharm, mpf
from itertools import combinations
mp.dps = 60
for r in results:
    k = r["k"]
    if k > 12:
        break
    ps = primes[:k]
    G = mpf(0)
    for rr in range(k + 1):
        for sub in combinations(ps, rr):
            d = 1
            for p in sub:
                d *= p
            G += (-1) ** rr * mharm(d)
    Gf = Fraction(mp.nstr(G, 55, strip_zeros=False))
    lo, hi = Fraction(r["G_lo"]), Fraction(r["G_hi"])
    dist = max(Fraction(0), lo - Gf, Gf - hi)
    assert dist < Fraction(1, 10 ** 50), (k, float(dist))
    r["bruteforce_60dps_distance"] = float(dist)
print("brute force over all 2^k divisors (60 digits) for k = 2..12: all within 1e-50 of my enclosures")
assert all(r["sign_G"] != 0 for r in results)
assert all(r["sign_G"] == (1 if r["k"] < first_neg else -1) for r in results)
assert all(r["inside_route_a"] and r["inside_route_b"] and r["sign_agrees_route_a"] and r["sign_agrees_route_b"] for r in results)
print("ALL %d enclosures sign-definite, contained in BOTH certificate enclosures, signs agree; G>0 for 2<=k<=%d, G<0 for %d<=k<=200"
      % (len(results), first_neg - 1, first_neg))
json.dump({"method": __doc__, "parameters": {"m": M, "D": D, "KMAX": KMAX, "mpfr_bits": PREC, "tail_exponent": S_TAIL},
           "gamma": [dec(GAMMA.to_frac()[0], 100, False), dec(GAMMA.to_frac()[1], 100, True)],
           "first_k_with_G_negative": first_neg, "results": results}, open(HERE + "/results_routeC.json", "w"), indent=1)
print("wrote results_routeC.json   [%.1fs]" % (time.time() - t0))
