#!/usr/bin/env python3
"""
Route (a): rigorous sign certificate for a(p_k#) (OEIS A067857) via mpmath.iv interval
arithmetic on the closed form.

Notation.  n = p_k# (product of the first k primes), k >= 2.
    F(n) = sum_{d|n} mu(n/d) H_d,      a(n) = n! F(n),
    G(k) := (-1)^k F(n) = sum_{d|n} (-1)^{omega(d)} H_d.
Conjecture (Israel): a(n) < 0 <=> omega(n) odd and >= 3.  At n = p_k#, k >= 3, it is
equivalent to G(k) > 0.  (G(k) < 0 with k even violates the "only if" direction, G(k) < 0
with k odd violates the "if" direction.)

Closed form (proved in README / JSON "justification"):
  eps(d) := H_d - ln d - gamma = A_m(d) + R_m(d),
  A_m(d)  = 1/(2d) - sum_{j=1}^m B_{2j} / (2j d^{2j}),
  |R_m(d)| < |B_{2m+2}| / ((2m+2) d^{2m+2})          [DLMF 5.11(ii); Alzer 1997, Thm 8]
  G(k) = 1/2 P_1 - sum_{j=1}^m B_{2j}/(2j) P_{2j}
         + sum_{d|n, d<=D} (-1)^{omega(d)} (eps(d) - A_m(d))
         + E,      |E| <= |B_{2m+2}|/(2m+2) * ( prod_{p|n}(1+p^{-(2m+2)}) - sum_{d|n, d<=D} d^{-(2m+2)} ),
  where P_s = prod_{p|n} (1 - p^{-s}).
Here m = 4 (s = 1, 2, 4, 6, 8; tail exponent 10), D = 10^5, 256-bit interval arithmetic.

Transcendental constants: ln p is enclosed by the atanh series with an explicit tail bound
(evaluated in interval arithmetic); gamma is enclosed from the exact interval H_D and the
same remainder lemma at x = D.  mpmath's iv.log / iv.euler are only used as consistency checks.
"""
import json, os, sys, time
from fractions import Fraction
from mpmath import iv, mp
from mpmath.libmp import to_rational
import sympy

HERE = os.path.dirname(os.path.abspath(__file__))
PREC = 256          # bits
M = 4               # Euler-Maclaurin terms: s = 1,2,4,6,8 ; tail exponent 2M+2 = 10
D = 10 ** 5         # small divisors d <= D treated exactly
KMAX = 200          # certify k = 2..KMAX
JSER = 90           # atanh-series terms for ln p  (t <= 1/3: tail <= 3^{-181}/181 < 1e-88)
OUT_DIGITS = 60     # decimal digits of enclosures written to the JSON

iv.prec = PREC
mp.prec = PREC
t0 = time.time()

# ---------------------------------------------------------------- helpers
def endpoints(x):
    """Exact rational endpoints of an iv.mpf (or degenerate value)."""
    a, b = x._mpi_
    return Fraction(*to_rational(a)), Fraction(*to_rational(b))

def fmt_fixed(q, digits):
    sign = '-' if q < 0 else ''
    s = str(abs(q)).rjust(digits + 1, '0')
    return sign + s[:-digits] + '.' + s[-digits:]

def dec_floor(fr, digits=OUT_DIGITS):
    return fmt_fixed((fr.numerator * 10 ** digits) // fr.denominator, digits)

def dec_ceil(fr, digits=OUT_DIGITS):
    return fmt_fixed(-((-fr.numerator * 10 ** digits) // fr.denominator), digits)

def enclosure_strings(x, digits=OUT_DIGITS):
    lo, hi = endpoints(x)
    return [dec_floor(lo, digits), dec_ceil(hi, digits)]

def ivfrac(fr):
    """iv.mpf enclosing an exact Fraction."""
    return iv.mpf(fr.numerator) / iv.mpf(fr.denominator)

def hull(lo, hi):
    return iv.mpf([lo, hi])

# ---------------------------------------------------------------- primes, Bernoulli numbers
primes = list(sympy.primerange(2, 10 ** 4))[:KMAX]
assert len(primes) == KMAX and primes[91] == 479 and primes[92] == 487
B = {2 * j: Fraction(int(sympy.bernoulli(2 * j).p), int(sympy.bernoulli(2 * j).q)) for j in range(1, M + 2)}
assert B[2] == Fraction(1, 6) and B[4] == Fraction(-1, 30) and B[10] == Fraction(5, 66)
S_TAIL = 2 * M + 2
C_TAIL = abs(B[S_TAIL]) / S_TAIL                    # |B_{2m+2}|/(2m+2) = 1/132 for m=4
coef = {1: Fraction(1, 2)}                          # G = sum_s coef[s] * P_s + corrections
for j in range(1, M + 1):
    coef[2 * j] = -B[2 * j] / (2 * j)
S_LIST = sorted(coef)

# ---------------------------------------------------------------- ln p via atanh series (interval)
def atanh_series_iv(t):
    """Rigorous enclosure of atanh(t) = sum_{j>=0} t^{2j+1}/(2j+1) for an interval t in [0, 1/3].
    Partial sum of JSER terms plus tail in [0, t^{2J+1} / ((2J+1)(1 - t^2))]."""
    t2 = t * t
    term = t
    s = iv.mpf(0)
    for j in range(JSER):
        s += term / (2 * j + 1)
        term = term * t2
    tail_hi = (term / ((2 * JSER + 1) * (1 - t2))).b      # term == t^{2J+1} here
    return s + hull(0, tail_hi)

LN2 = 2 * atanh_series_iv(iv.mpf(1) / iv.mpf(3))            # ln 2 = 2 atanh(1/3)

def ln_int_iv(p):
    """Rigorous enclosure of ln p for an integer p >= 1: p = 2^e y, y in [1,2), ln y = 2 atanh((y-1)/(y+1))."""
    if p == 1:
        return iv.mpf(0)
    e = p.bit_length() - 1
    t = iv.mpf(p - 2 ** e) / iv.mpf(p + 2 ** e)             # in [0, 1/3)
    return e * LN2 + 2 * atanh_series_iv(t)

LN = {p: ln_int_iv(p) for p in primes}
# consistency check against mpmath's own log (independent algorithm); not used in the certificate
for p in primes:
    lo, hi = endpoints(LN[p])
    mlo, mhi = endpoints(iv.log(iv.mpf(p)))
    assert lo <= mhi and mlo <= hi, ("ln %d series enclosure inconsistent with iv.log" % p)
    assert (hi - lo) < Fraction(1, 10 ** 70)
print("ln p enclosures (series, interval arithmetic) computed and consistent with iv.log for %d primes; width(ln 2) = %.3e"
      % (len(primes), float(endpoints(LN[2])[1] - endpoints(LN[2])[0])))

# ---------------------------------------------------------------- squarefree divisors d <= D of p_KMAX#
# entry: (d, omega, gpf_index (1-based, 0 for d=1), tuple of prime indices)
divs = []
def dfs(start, d, idxs):
    divs.append((d, len(idxs), (idxs[-1] + 1) if idxs else 0, tuple(idxs)))
    for i in range(start, KMAX):
        nd = d * primes[i]
        if nd > D:
            break
        dfs(i + 1, nd, idxs + [i])
dfs(0, 1, [])
divs.sort()
div_set = {d: (w, j, idxs) for d, w, j, idxs in divs}
print("squarefree p_%d-smooth d <= %d: %d divisors" % (KMAX, D, len(divs)))

# ---------------------------------------------------------------- interval harmonic numbers H_d, d <= D
Hd = {}
H = iv.mpf(0)
one = iv.mpf(1)
for i in range(1, D + 1):
    H = H + one / iv.mpf(i)
    if i in div_set or i == D:
        Hd[i] = H
HD = Hd[D]
print("interval H_d for all divisors and H_D computed; width(H_D) = %.3e   [%.1fs]"
      % (float(endpoints(HD)[1] - endpoints(HD)[0]), time.time() - t0))

# ---------------------------------------------------------------- A_m(d) and the tail power d^{-(2m+2)}
def A_m_iv(d):
    dd = iv.mpf(d)
    inv2 = one / (dd * dd)
    powers = {2: inv2}
    powers[4] = inv2 * inv2
    powers[6] = powers[4] * inv2
    powers[8] = powers[4] * powers[4]
    val = one / (2 * dd)
    for j in range(1, M + 1):
        val -= ivfrac(B[2 * j] / (2 * j)) * powers[2 * j]
    return val, powers[8] * inv2                              # (A_m(d), d^{-10})

# ---------------------------------------------------------------- gamma from H_D (remainder lemma at x = D)
LN_D = 5 * (LN[2] + LN[5])                                     # ln 10^5
A_D, _ = A_m_iv(D)
GAMMA = HD - LN_D - A_D + hull(-ivfrac(C_TAIL / Fraction(D) ** S_TAIL), ivfrac(C_TAIL / Fraction(D) ** S_TAIL))
glo, ghi = endpoints(GAMMA)
elo, ehi = endpoints(+iv.euler)
assert glo <= ehi and elo <= ghi, "derived gamma enclosure inconsistent with iv.euler"
print("gamma enclosure derived from H_D: width %.3e, consistent with iv.euler" % float(ghi - glo))

# ---------------------------------------------------------------- per-divisor corrections, bucketed by gpf index
bucketV = [iv.mpf(0) for _ in range(KMAX + 1)]     # sum of (-1)^omega (eps(d) - A_m(d)) over d with gpf index j
bucketS = [iv.mpf(0) for _ in range(KMAX + 1)]     # sum of d^{-10}
bucket_count = [0] * (KMAX + 1)
for d, w, j, idxs in divs:
    ln_d = iv.mpf(0)
    for i in idxs:
        ln_d += LN[primes[i]]
    A, dm10 = A_m_iv(d)
    v = Hd[d] - ln_d - GAMMA - A                              # eps(d) - A_m(d) = R_m(d)
    if w % 2:
        v = -v
    bucketV[j] += v
    bucketS[j] += dm10
    bucket_count[j] += 1
print("per-divisor corrections done   [%.1fs]" % (time.time() - t0))

# ---------------------------------------------------------------- per-k certificates
certs = []
P = {s: one for s in S_LIST}      # running Euler products prod_{i<=k} (1 - p_i^{-s})
Q = one                            # running prod (1 + p_i^{-10})
corr = iv.mpf(0)
tail_small = iv.mpf(0)
ndiv = 0
corr += bucketV[0]; tail_small += bucketS[0]; ndiv += bucket_count[0]    # d = 1
sign_changes = []
first_negative = None
for k in range(1, KMAX + 1):
    p = primes[k - 1]
    pp = iv.mpf(p)
    for s in S_LIST:
        P[s] = P[s] * (one - one / pp ** s)
    Q = Q * (one + one / pp ** S_TAIL)
    corr += bucketV[k]; tail_small += bucketS[k]; ndiv += bucket_count[k]
    if k < 2:
        continue
    closed = iv.mpf(0)
    for s in S_LIST:
        closed += ivfrac(coef[s]) * P[s]
    tail_rest_hi = max(Fraction(0), endpoints(Q - tail_small)[1])       # >= sum_{d|n, d>D} d^{-10}
    T = C_TAIL * tail_rest_hi                                            # exact rational bound on |E|
    G = closed + corr + hull(-ivfrac(T), ivfrac(T))
    glo, ghi = endpoints(G)
    if ghi < 0:
        sgnG = -1
    elif glo > 0:
        sgnG = 1
    else:
        sgnG = 0
    a_sign = (-1) ** k * sgnG            # sign of F(n) = sign of a(n)
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
        "euler_products": {str(s): enclosure_strings(P[s]) for s in S_LIST},
        "closed_form": enclosure_strings(closed),
        "small_divisor_correction": enclosure_strings(corr),
        "tail_divisor_sum_upper_bound": dec_ceil(tail_rest_hi, 80),
        "tail_bound_T": dec_ceil(T, 80),
        "G_enclosure": enclosure_strings(G),
        "G_enclosure_float": [float(glo), float(ghi)],
        "sign_G": sgnG,
        "sign_a": a_sign,
        "conjecture_predicts_a_negative": predicted_negative,
        "status": status,
    })
    if k >= 85 and k <= 100 or k <= 5 or k % 25 == 0:
        print("k=%3d p_k=%4d  G in [%+.6e, %+.6e]  sign %+d  %s" % (k, p, float(glo), float(ghi), sgnG, status))

print("first k with certified G(k) < 0:", first_negative)
assert all(c["sign_G"] != 0 for c in certs), "some enclosure is not sign-definite"
assert all(c["sign_G"] == 1 for c in certs if c["k"] <= first_negative - 1)
assert all(c["sign_G"] == -1 for c in certs if c["k"] >= first_negative)

justification = [
    "L1 (digamma remainder; DLMF 5.11(ii), Alzer 1997 Math. Comp. 66, Thm 8): for real x>0 and integer m>=0, "
    "psi(x) = ln x - 1/(2x) - sum_{j=1}^m B_{2j}/(2j x^{2j}) + R_m(x) with R_m(x) of the sign of -B_{2m+2} and "
    "|R_m(x)| < |B_{2m+2}|/((2m+2) x^{2m+2}).  Since H_d = psi(d) + 1/d + gamma, eps(d) := H_d - ln d - gamma = A_m(d) + R_m(d).",
    "L2: for squarefree n with k = omega(n) >= 2: sum_{d|n} (-1)^{omega(d)} = 0 and sum_{d|n} (-1)^{omega(d)} ln d = "
    "sum_{p|n} ln p * sum_{d: p|d|n} (-1)^{omega(d)} = -(1-1)^{k-1} sum_p ln p = 0; hence G(k) = sum_{d|n} (-1)^{omega(d)} eps(d).",
    "L3: sum_{d|n} (-1)^{omega(d)} d^{-s} = prod_{p|n} (1 - p^{-s}) and sum_{d|n} d^{-s} = prod_{p|n} (1 + p^{-s}) (multiplicativity).",
    "Combining: G(k) = 1/2 P_1 - sum_j B_{2j}/(2j) P_{2j} + sum_{d<=D} (-1)^{omega(d)} (eps(d) - A_m(d)) + E, "
    "|E| <= sum_{d|n, d>D} |R_m(d)| <= |B_{2m+2}|/(2m+2) * (prod_{p|n}(1+p^{-(2m+2)}) - sum_{d|n,d<=D} d^{-(2m+2)}).",
    "ln p: p = 2^e y, y in [1,2), ln y = 2 atanh(t), t = (y-1)/(y+1) in [0,1/3]; atanh(t) = sum_{j>=0} t^{2j+1}/(2j+1); "
    "the tail after J terms lies in [0, t^{2J+1}/((2J+1)(1-t^2))]; ln 2 = 2 atanh(1/3).",
    "gamma = H_D - ln D - eps(D) = H_D - ln D - A_m(D) - R_m(D), with |R_m(D)| < |B_{2m+2}|/((2m+2) D^{2m+2}) by L1.",
    "Arithmetic: mpmath.iv with directed rounding (basic operations are correctly rounded); H_d as running interval sums; "
    "all interval endpoints are converted exactly to rationals and rounded outward to decimal strings.",
    "F(n) = (-1)^k G(k) and a(n) = n! F(n), n! > 0; so sign a(p_k#) = (-1)^k sign G(k).",
]
out = {
    "sequence": "A067857",
    "route": "a: mpmath.iv interval arithmetic on the closed form",
    "parameters": {"precision_bits": PREC, "m": M, "s_list": S_LIST, "tail_exponent": S_TAIL, "D": D, "KMAX": KMAX,
                   "atanh_series_terms": JSER, "output_decimal_digits": OUT_DIGITS},
    "primes": primes,
    "bernoulli": {str(s): "%d/%d" % (B[s].numerator, B[s].denominator) for s in sorted(B)},
    "coefficients": {str(s): "%d/%d" % (coef[s].numerator, coef[s].denominator) for s in S_LIST},
    "tail_constant": "%d/%d" % (C_TAIL.numerator, C_TAIL.denominator),
    "constants": {"gamma": enclosure_strings(GAMMA, 80), "ln_p": {str(p): enclosure_strings(LN[p], 80) for p in primes},
                  "ln2": enclosure_strings(LN2, 80)},
    "H_D": enclosure_strings(HD, 80),
    "num_small_divisors_total": len(divs),
    "first_k_with_G_negative": first_negative,
    "justification": justification,
    "certificates": certs,
}
with open(os.path.join(HERE, "cert_route_a.json"), "w") as f:
    json.dump(out, f, indent=1)
print("wrote cert_route_a.json  [%.1fs]" % (time.time() - t0))
