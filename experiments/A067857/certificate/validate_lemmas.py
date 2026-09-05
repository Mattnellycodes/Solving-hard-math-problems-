#!/usr/bin/env python3
"""
Numerical sanity checks (NOT proofs) of the lemmas used by the certificates, at 80-100 digits.

  L1  (digamma remainder, DLMF 5.11(ii) / Alzer 1997 Thm 8):
      psi(x) = ln x - 1/(2x) - sum_{j=1}^m B_{2j}/(2j x^{2j}) + R_m(x),
      sign R_m(x) = sign(-B_{2m+2}),  |R_m(x)| < |B_{2m+2}| / ((2m+2) x^{2m+2}),  for real x > 0.
      Also H_d = psi(d) + 1/d + gamma.
  L2  For squarefree n with omega(n) >= 2: sum_{d|n} (-1)^{omega(d)} = 0 and sum_{d|n} (-1)^{omega(d)} ln d = 0.
  L3  sum_{d|n} (-1)^{omega(d)} d^{-s} = prod_{p|n} (1 - p^{-s}),  sum_{d|n} d^{-s} = prod_{p|n} (1 + p^{-s})  (exact check).
  C   The rational atanh-series enclosures of ln p and the H_D-derived gamma enclosure contain
      mpmath's 100-digit values.
"""
import sys
from fractions import Fraction
from math import comb
from mpmath import mp, mpf, psi, log, harmonic, euler, bernoulli as mpbern

sys.set_int_max_str_digits(0)
mp.dps = 100

def bernoulli_even(nmax):
    Bs = [Fraction(1)]
    for n in range(1, nmax + 1):
        Bs.append(-sum(comb(n + 1, j) * Bs[j] for j in range(n)) / (n + 1))
    return {2 * j: Bs[2 * j] for j in range(1, nmax // 2 + 1)}
B = bernoulli_even(20)
for s in B:
    assert abs(mpbern(s) - mpf(B[s].numerator) / B[s].denominator) < mpf(10) ** -90
print("Bernoulli numbers from the recurrence agree with mpmath:", {s: str(B[s]) for s in sorted(B) if s <= 14})

# ---- L1
bad = 0
tests = 0
for x in [1, 2, 3, 4, 5, 7, 10, 30, 100, 1000, 10 ** 4, 10 ** 5, mpf('0.5'), mpf('1.5'), mpf('2.25')]:
    x = mpf(x)
    for m in range(0, 9):
        S = log(x) - 1 / (2 * x) - sum(mpf(B[2 * j].numerator) / B[2 * j].denominator / (2 * j * x ** (2 * j)) for j in range(1, m + 1))
        R = psi(0, x) - S
        bound = abs(mpf(B[2 * m + 2].numerator) / B[2 * m + 2].denominator) / ((2 * m + 2) * x ** (2 * m + 2))
        sgn_expected = -1 if B[2 * m + 2] > 0 else 1
        tests += 1
        if not (R * sgn_expected > 0 and abs(R) < bound):
            bad += 1
            print("   L1 FAILS at x=%s m=%d: R=%s bound=%s" % (x, m, R, bound))
print("L1 (sign and magnitude of the digamma remainder): %d/%d cases pass" % (tests - bad, tests))
for d in [1, 2, 3, 10, 1000, 10 ** 5]:
    assert abs(harmonic(d) - (psi(0, mpf(d)) + mpf(1) / d + euler)) < mpf(10) ** -95
print("H_d = psi(d) + 1/d + gamma verified numerically for several d")

# ---- L2, L3
def first_primes(k):
    ps, q = [], 2
    while len(ps) < k:
        if all(q % p for p in ps):
            ps.append(q)
        q += 1
    return ps
def divisors_with_omega(ps):
    divs = [(1, 0)]
    for p in ps:
        divs += [(d * p, w + 1) for d, w in divs]
    return divs
for k in range(2, 11):
    ps = first_primes(k)
    divs = divisors_with_omega(ps)
    s0 = sum((-1) ** w for d, w in divs)
    s1 = sum(((-1) ** w) * log(mpf(d)) for d, w in divs)
    assert s0 == 0 and abs(s1) < mpf(10) ** -90, (k, s0, s1)
    for s in (1, 2, 4, 6, 8, 10, 12, 14):
        lhs = sum(Fraction((-1) ** w, d ** s) for d, w in divs)
        rhs = Fraction(1)
        for p in ps:
            rhs *= 1 - Fraction(1, p ** s)
        assert lhs == rhs
        lhs2 = sum(Fraction(1, d ** s) for d, w in divs)
        rhs2 = Fraction(1)
        for p in ps:
            rhs2 *= 1 + Fraction(1, p ** s)
        assert lhs2 == rhs2
print("L2 (sum (-1)^omega = 0 and sum (-1)^omega ln d = 0) checked numerically for k=2..10;")
print("L3 (Euler-product identities) checked EXACTLY for k=2..10, s in {1,2,...,14}")
# L2 proof sketch check for k=1 (where it fails): n = 2: sum (-1)^omega ln d = -ln 2 != 0
assert abs(sum(((-1) ** w) * log(mpf(d)) for d, w in divisors_with_omega([2])) + log(2)) < mpf(10) ** -90
print("   (and for k = 1 the ln-sum is -ln p, as expected: the closed form needs k >= 2)")

# ---- C: constants
JSER = 100
def atanh_enclosure(t):
    t2, term, s = t * t, t, Fraction(0)
    for j in range(JSER):
        s += term / (2 * j + 1)
        term *= t2
    return s, s + term / ((2 * JSER + 1) * (1 - t2))
LN2 = tuple(2 * x for x in atanh_enclosure(Fraction(1, 3)))
def ln_int_enclosure(p):
    e = p.bit_length() - 1
    lo, hi = atanh_enclosure(Fraction(p - 2 ** e, p + 2 ** e))
    return (e * LN2[0] + 2 * lo, e * LN2[1] + 2 * hi)
def contains(enc, val):
    lo, hi = enc
    return mpf(lo.numerator) / lo.denominator <= val <= mpf(hi.numerator) / hi.denominator
worst = 0
for p in first_primes(200) + [10, 1000, 4096, 4097]:
    enc = ln_int_enclosure(p)
    assert contains(enc, log(mpf(p))), p
    worst = max(worst, float(enc[1] - enc[0]))
print("C: rational atanh-series enclosures of ln p contain mpmath's 100-digit log(p) for the first 200 primes (max width %.1e)" % worst)

# gamma from H_D with the remainder lemma, D = 10^4, m = 6 and D = 10^5, m = 4
for (D, m) in [(10 ** 4, 6), (10 ** 5, 4)]:
    N, L = 0, 1
    spf = list(range(D + 1))
    for i in range(2, int(D ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, D + 1, i):
                if spf[j] == j:
                    spf[j] = i
    for d in range(1, D + 1):
        if d > 1:
            q, x = spf[d], d
            while x % q == 0:
                x //= q
            if x == 1:
                L *= q; N *= q
        N += L // d
    HD = Fraction(N, L)
    A = 1 / (2 * Fraction(D)) - sum(B[2 * j] / (2 * j * Fraction(D) ** (2 * j)) for j in range(1, m + 1))
    C = abs(B[2 * m + 2]) / (2 * m + 2) / Fraction(D) ** (2 * m + 2)
    a10 = len(str(D)) - 1
    lnD = (a10 * (LN2[0] + ln_int_enclosure(5)[0]), a10 * (LN2[1] + ln_int_enclosure(5)[1]))
    g = (HD - A - lnD[1] - C, HD - A - lnD[0] + C)
    assert contains(g, euler)
    print("C: gamma enclosure from H_%d (m=%d) contains mpmath's 100-digit euler constant; width %.1e" % (D, m, float(g[1] - g[0])))
print("all lemma sanity checks passed")
