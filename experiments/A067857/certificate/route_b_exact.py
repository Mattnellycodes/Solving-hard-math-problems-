#!/usr/bin/env python3
"""
Route (b): rigorous sign certificate for a(p_k#) (OEIS A067857) using EXACT RATIONAL
arithmetic only (python fractions / integers; no mpmath, no sympy).

Same closed form as route (a) (see docstring there and the JSON "justification"), but with
m = 6 (s = 1,2,4,6,8,10,12; tail exponent 14) and D = 10^4:

  G(k) = (-1)^k F(p_k#) = Q(k) - [ sum_p c_p(k) ln p + c_0(k) gamma ] + E(k),

where everything in Q(k) is an exact rational:
  Q(k) = 1/2 P_1 - sum_{j<=m} B_{2j}/(2j) P_{2j}
         + sum_{d | p_k#, d <= D} (-1)^{omega(d)} ( H_d - A_m(d) ),
  c_p(k) = sum_{d | p_k#, d <= D, p | d} (-1)^{omega(d)},     c_0(k) = sum_{d | p_k#, d <= D} (-1)^{omega(d)},
  |E(k)| <= T(k) := |B_{2m+2}|/(2m+2) * ( prod_{i<=k} (1 + p_i^{-(2m+2)}) - sum_{d | p_k#, d<=D} d^{-(2m+2)} )
(this is exactly sum_{d<=D} (-1)^omega (eps(d) - A_m(d)) with eps(d) = H_d - ln d - gamma,
 the transcendental parts collected into the integer coefficients c_p, c_0).
Only ln p and gamma are enclosed by intervals (rational endpoints):
  ln p from the atanh series with explicit tail bound, gamma from the exact H_D and the
  remainder lemma L1 at x = D.
"""
import hashlib, json, os, sys, time
from fractions import Fraction
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.set_int_max_str_digits(0)      # exact rationals here have > 4300-digit numerators
M = 6
D = 10 ** 4
KMAX = 200
JSER = 100
OUT_DIGITS = 80
t0 = time.time()

# ---------------------------------------------------------------- helpers
def fmt_fixed(q, digits):
    sign = '-' if q < 0 else ''
    s = str(abs(q)).rjust(digits + 1, '0')
    return sign + s[:-digits] + '.' + s[-digits:]

def dec_floor(fr, digits=OUT_DIGITS):
    return fmt_fixed((fr.numerator * 10 ** digits) // fr.denominator, digits)

def dec_ceil(fr, digits=OUT_DIGITS):
    return fmt_fixed(-((-fr.numerator * 10 ** digits) // fr.denominator), digits)

def enc(lo, hi, digits=OUT_DIGITS):
    assert lo <= hi
    return [dec_floor(lo, digits), dec_ceil(hi, digits)]

def frac_hash(fr):
    return hashlib.sha256(("%d/%d" % (fr.numerator, fr.denominator)).encode()).hexdigest()[:16]

# ---------------------------------------------------------------- primes (own sieve), Bernoulli (recurrence)
LIM = 2000
sieve = bytearray([1]) * (LIM + 1)
sieve[0] = sieve[1] = 0
for i in range(2, int(LIM ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
primes = [i for i in range(LIM + 1) if sieve[i]][:KMAX]
assert len(primes) == KMAX and primes[0] == 2 and primes[91] == 479 and primes[92] == 487 and primes[199] == 1223

def bernoulli_list(nmax):
    Bs = [Fraction(1)]
    for n in range(1, nmax + 1):
        Bs.append(-sum(comb(n + 1, j) * Bs[j] for j in range(n)) / (n + 1))
    return Bs
Bl = bernoulli_list(2 * M + 2)
B = {2 * j: Bl[2 * j] for j in range(1, M + 2)}
assert B[2] == Fraction(1, 6) and B[4] == Fraction(-1, 30) and B[12] == Fraction(-691, 2730) and B[14] == Fraction(7, 6)
S_TAIL = 2 * M + 2
C_TAIL = abs(B[S_TAIL]) / S_TAIL
coef = {1: Fraction(1, 2)}
for j in range(1, M + 1):
    coef[2 * j] = -B[2 * j] / (2 * j)
S_LIST = sorted(coef)

def A_m(d):
    d = Fraction(d)
    return 1 / (2 * d) - sum(B[2 * j] / (2 * j * d ** (2 * j)) for j in range(1, M + 1))

# ---------------------------------------------------------------- ln p enclosures (exact rational series)
def atanh_enclosure(t):
    """[lo, hi] rational enclosure of atanh(t) for rational 0 <= t <= 1/3."""
    assert 0 <= t <= Fraction(1, 3)
    t2 = t * t
    term = t
    s = Fraction(0)
    for j in range(JSER):
        s += term / (2 * j + 1)
        term *= t2
    tail_hi = term / ((2 * JSER + 1) * (1 - t2))
    return s, s + tail_hi

LN2 = tuple(2 * x for x in atanh_enclosure(Fraction(1, 3)))

def ln_int_enclosure(p):
    if p == 1:
        return (Fraction(0), Fraction(0))
    e = p.bit_length() - 1
    lo, hi = atanh_enclosure(Fraction(p - 2 ** e, p + 2 ** e))
    return (e * LN2[0] + 2 * lo, e * LN2[1] + 2 * hi)

LN = {p: ln_int_enclosure(p) for p in primes}
print("ln p enclosures: width(ln 2) = %.3e, width(ln 1223) = %.3e" % (float(LN[2][1] - LN[2][0]), float(LN[1223][1] - LN[1223][0])))

# ---------------------------------------------------------------- divisors d <= D of p_KMAX#
divs = []
def dfs(start, d, idxs):
    divs.append((d, len(idxs), (idxs[-1] + 1) if idxs else 0, tuple(idxs)))
    for i in range(start, KMAX):
        nd = d * primes[i]
        if nd > D:
            break
        dfs(i + 1, nd, idxs + [i])
dfs(0, 1, [])
div_info = {d: (w, j, idxs) for d, w, j, idxs in divs}
print("squarefree p_%d-smooth d <= %d: %d divisors" % (KMAX, D, len(divs)))

# ---------------------------------------------------------------- exact H_d via lcm accumulator, bucketed sums
spf = list(range(D + 1))
for i in range(2, int(D ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, D + 1, i):
            if spf[j] == j:
                spf[j] = i

def prime_power_base(d):
    """Return q if d = q^a for a prime q (a >= 1), else None."""
    q = spf[d]
    x = d
    while x % q == 0:
        x //= q
    return q if x == 1 else None

bucketH = [Fraction(0) for _ in range(KMAX + 1)]      # sum (-1)^omega H_d
bucketA = [Fraction(0) for _ in range(KMAX + 1)]      # sum (-1)^omega A_m(d)
bucketS = [Fraction(0) for _ in range(KMAX + 1)]      # sum d^{-(2m+2)}
bucketC0 = [0] * (KMAX + 1)                            # sum (-1)^omega
bucketCp = [dict() for _ in range(KMAX + 1)]           # prime index -> sum (-1)^omega over d containing it
bucket_count = [0] * (KMAX + 1)
N, L = 0, 1                                            # H_d = N / L, L = lcm(1..d)
HD = None
for d in range(1, D + 1):
    if d > 1:
        q = prime_power_base(d)
        if q is not None:
            L *= q
            N *= q
    N += L // d
    if d in div_info:
        w, j, idxs = div_info[d]
        sg = -1 if w % 2 else 1
        Hd = Fraction(N, L)
        bucketH[j] += sg * Hd
        bucketA[j] += sg * A_m(d)
        bucketS[j] += Fraction(1, d ** S_TAIL)
        bucketC0[j] += sg
        for i in idxs:
            bucketCp[j][i] = bucketCp[j].get(i, 0) + sg
        bucket_count[j] += 1
    if d == D:
        HD = Fraction(N, L)
assert HD is not None
print("exact H_d for all divisors and H_D done (denominator of H_D has %d digits)   [%.1fs]"
      % (len(str(HD.denominator)), time.time() - t0))

# ---------------------------------------------------------------- gamma from H_D
LN_D = (4 * (LN[2][0] + LN[5][0]), 4 * (LN[2][1] + LN[5][1]))           # ln 10^4
RD = C_TAIL / Fraction(D) ** S_TAIL
GAMMA = (HD - A_m(D) - LN_D[1] - RD, HD - A_m(D) - LN_D[0] + RD)
print("gamma in [%s, %s] (width %.3e)" % (dec_floor(GAMMA[0], 40), dec_ceil(GAMMA[1], 40), float(GAMMA[1] - GAMMA[0])))
assert Fraction("0.5772156649015328606065120900824024310421") < GAMMA[0] and GAMMA[1] < Fraction("0.5772156649015328606065120900824024310422")

# ---------------------------------------------------------------- per-k certificates
certs = []
P = {s: Fraction(1) for s in S_LIST}
Qt = Fraction(1)
Hsum = Fraction(0); Asum = Fraction(0); Ssum = Fraction(0); c0 = 0
cp = [0] * KMAX
ndiv = 0
def absorb(j):
    global Hsum, Asum, Ssum, c0, ndiv
    Hsum += bucketH[j]; Asum += bucketA[j]; Ssum += bucketS[j]; c0 += bucketC0[j]; ndiv += bucket_count[j]
    for i, v in bucketCp[j].items():
        cp[i] += v
absorb(0)
first_negative = None
for k in range(1, KMAX + 1):
    p = primes[k - 1]
    for s in S_LIST:
        P[s] *= (1 - Fraction(1, p ** s))
    Qt *= (1 + Fraction(1, p ** S_TAIL))
    absorb(k)
    if k < 2:
        continue
    closed = sum(coef[s] * P[s] for s in S_LIST)
    Qk = closed + Hsum - Asum                                     # exact rational
    # transcendental part  X = sum_p c_p ln p + c0 gamma  (interval)
    Xlo = Xhi = Fraction(0)
    for i in range(k):
        c = cp[i]
        if c > 0:
            Xlo += c * LN[primes[i]][0]; Xhi += c * LN[primes[i]][1]
        elif c < 0:
            Xlo += c * LN[primes[i]][1]; Xhi += c * LN[primes[i]][0]
    if c0 > 0:
        Xlo += c0 * GAMMA[0]; Xhi += c0 * GAMMA[1]
    elif c0 < 0:
        Xlo += c0 * GAMMA[1]; Xhi += c0 * GAMMA[0]
    tail_rest = Qt - Ssum                                         # = sum_{d | p_k#, d > D} d^{-14}  (exact)
    assert tail_rest >= 0
    T = C_TAIL * tail_rest
    Glo = Qk - Xhi - T
    Ghi = Qk - Xlo + T
    sgnG = -1 if Ghi < 0 else (1 if Glo > 0 else 0)
    a_sign = (-1) ** k * sgnG
    predicted_negative = (k % 2 == 1 and k >= 3)
    if sgnG == 0:
        status = "UNDECIDED"
    elif (a_sign < 0) == predicted_negative:
        status = "conjecture holds at n=p_k#"
    else:
        status = ("VIOLATION of 'only if' (a(n)<0 with omega even)" if k % 2 == 0
                  else "VIOLATION of 'if' (a(n)>0 with omega odd >= 3)")
    if sgnG < 0 and first_negative is None:
        first_negative = k
    certs.append({
        "k": k, "p_k": p, "omega": k, "omega_parity": "even" if k % 2 == 0 else "odd",
        "num_small_divisors_used": ndiv,
        "euler_products_exact": {str(s): enc(P[s], P[s]) for s in S_LIST},
        "closed_form_exact": enc(closed, closed),
        "Hsum_exact": enc(Hsum, Hsum),
        "Asum_exact": enc(Asum, Asum),
        "Q_exact": enc(Qk, Qk),
        "Q_exact_sha256_16": frac_hash(Qk),
        "c0": c0,
        "c_p": {str(primes[i]): cp[i] for i in range(k) if cp[i] != 0},
        "transcendental_part": enc(Xlo, Xhi),
        "tail_divisor_sum_exact": enc(tail_rest, tail_rest),
        "tail_bound_T": dec_ceil(T),
        "G_enclosure": enc(Glo, Ghi),
        "G_enclosure_float": [float(Glo), float(Ghi)],
        "sign_G": sgnG, "sign_a": a_sign,
        "conjecture_predicts_a_negative": predicted_negative,
        "status": status,
    })
    if 85 <= k <= 100 or k <= 5 or k % 25 == 0:
        print("k=%3d p_k=%4d  G in [%+.6e, %+.6e]  width %.1e  sign %+d  %s" % (k, p, float(Glo), float(Ghi), float(Ghi - Glo), sgnG, status))

print("first k with certified G(k) < 0:", first_negative, "  [%.1fs]" % (time.time() - t0))
assert all(c["sign_G"] != 0 for c in certs)
assert all(c["sign_G"] == 1 for c in certs if c["k"] < first_negative)
assert all(c["sign_G"] == -1 for c in certs if c["k"] >= first_negative)

justification = [
    "L1 (digamma remainder; DLMF 5.11(ii), Alzer 1997 Math. Comp. 66, Thm 8): for real x>0 and integer m>=0, "
    "psi(x) = ln x - 1/(2x) - sum_{j=1}^m B_{2j}/(2j x^{2j}) + R_m(x) with R_m(x) of the sign of -B_{2m+2} and "
    "|R_m(x)| < |B_{2m+2}|/((2m+2) x^{2m+2}).  Since H_d = psi(d) + 1/d + gamma, eps(d) := H_d - ln d - gamma = A_m(d) + R_m(d).",
    "L2: for squarefree n with k = omega(n) >= 2: sum_{d|n} (-1)^{omega(d)} = 0 and sum_{d|n} (-1)^{omega(d)} ln d = 0; "
    "hence G(k) = sum_{d|n} (-1)^{omega(d)} eps(d).",
    "L3: sum_{d|n} (-1)^{omega(d)} d^{-s} = prod_{p|n} (1 - p^{-s}) and sum_{d|n} d^{-s} = prod_{p|n} (1 + p^{-s}).",
    "G(k) = Q(k) - X(k) + E(k): Q(k) exact rational (closed form + sum_{d<=D} (-1)^omega (H_d - A_m(d))), "
    "X(k) = sum_p c_p ln p + c_0 gamma with integer c_p, c_0 (from the ln d + gamma parts of eps(d), d<=D), "
    "|E(k)| <= T(k) = |B_{2m+2}|/(2m+2) * (prod(1+p^{-(2m+2)}) - sum_{d<=D} d^{-(2m+2)}) by L1 and L3.",
    "ln p: p = 2^e y, y in [1,2), ln y = 2 atanh((y-1)/(y+1)), atanh(t) = sum t^{2j+1}/(2j+1), tail after J terms in "
    "[0, t^{2J+1}/((2J+1)(1-t^2))] for 0<=t<=1/3; ln 2 = 2 atanh(1/3).",
    "gamma = H_D - ln D - A_m(D) - R_m(D) with |R_m(D)| < |B_{2m+2}|/((2m+2) D^{2m+2}) (L1 at x = D = 10^4).",
    "All other arithmetic is exact (python integers / fractions); decimal strings are outward-rounded.",
    "sign a(p_k#) = sign F(p_k#) = (-1)^k sign G(k).",
]
out = {
    "sequence": "A067857",
    "route": "b: exact rational arithmetic (fractions), no floating point",
    "parameters": {"m": M, "s_list": S_LIST, "tail_exponent": S_TAIL, "D": D, "KMAX": KMAX,
                   "atanh_series_terms": JSER, "output_decimal_digits": OUT_DIGITS},
    "primes": primes,
    "bernoulli": {str(s): "%d/%d" % (B[s].numerator, B[s].denominator) for s in sorted(B)},
    "coefficients": {str(s): "%d/%d" % (coef[s].numerator, coef[s].denominator) for s in S_LIST},
    "tail_constant": "%d/%d" % (C_TAIL.numerator, C_TAIL.denominator),
    "constants": {"gamma": enc(*GAMMA), "ln_p": {str(p): enc(*LN[p]) for p in primes}, "ln2": enc(*LN2)},
    "H_D": enc(HD, HD),
    "H_D_sha256_16": frac_hash(HD),
    "num_small_divisors_total": len(divs),
    "first_k_with_G_negative": first_negative,
    "justification": justification,
    "certificates": certs,
}
with open(os.path.join(HERE, "cert_route_b.json"), "w") as f:
    json.dump(out, f, indent=1)
print("wrote cert_route_b.json  [%.1fs]" % (time.time() - t0))
