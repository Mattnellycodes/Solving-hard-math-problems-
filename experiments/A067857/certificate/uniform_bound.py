#!/usr/bin/env python3
"""
Bonus (rigorous, exact rationals): a uniform-in-k upper bound proving G(k) < 0 for ALL k >= k0.

With D = 1 the closed form (k >= 2) reads
    G(k) = sum_s c_s P_s(k) + (eps(1) - A_m(1)) + E(k),   |E(k)| <= C_m ( prod_{i<=k}(1+p_i^{-(2m+2)}) - 1 ),
    c_1 = 1/2, c_{2j} = -B_{2j}/(2j),  C_m = |B_{2m+2}|/(2m+2),  eps(1) = 1 - gamma,
    P_s(k) = prod_{i<=k} (1 - p_i^{-s}).
Facts (s >= 2):  P_s(k) is strictly decreasing in k;  P_s(k) >= prod_{all p}(1-p^{-s}) = 1/zeta(s);
                 prod_{i<=k}(1+p_i^{-s}) <= prod_{all p}(1+p^{-s}) = zeta(s)/zeta(2s).
Hence, with S+ = {s : c_s > 0} and S- = {s : c_s < 0},
    G(k) <= UB(k) := sum_{s in S+} c_s P_s(k) + sum_{s in S-} c_s / zeta(s) + (1 - gamma - A_m(1)) + C_m (zeta(2m+2)/zeta(4m+4) - 1),
and UB(k) is strictly decreasing in k.  So UB(k0) < 0 implies G(k) < 0 for every k >= k0, i.e.
    sign a(p_k#) = (-1)^{k+1} for all k >= k0
(a(p_k#) < 0 for even k >= k0 and a(p_k#) > 0 for odd k >= k0: infinitely many counterexamples to
both directions of the conjecture).  No Mertens-type estimate is used.

zeta(2n) = |B_{2n}| (2 pi)^{2n} / (2 (2n)!) with pi enclosed by Machin's formula
pi = 16 arctan(1/5) - 4 arctan(1/239) (alternating series: the sum lies between consecutive partial sums).
gamma: enclosure derived from H_D, D = 10^4, m = 6 (lemma L1), as in route (b).
"""
import json, os, sys
from fractions import Fraction
from math import comb, factorial

sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__))

def bernoulli_even(nmax):
    Bs = [Fraction(1)]
    for n in range(1, nmax + 1):
        Bs.append(-sum(comb(n + 1, j) * Bs[j] for j in range(n)) / (n + 1))
    return {2 * j: Bs[2 * j] for j in range(1, nmax // 2 + 1)}

B = bernoulli_even(32)

def arctan_enclosure(x, J=70):
    """x rational in (0,1]: arctan x lies between the partial sums S_J and S_{J+1} (alternating, decreasing terms)."""
    s = Fraction(0)
    term = x
    x2 = x * x
    S = []
    for j in range(J + 2):
        s += (-1) ** j * term / (2 * j + 1)
        S.append(s)
        term *= x2
    return (min(S[J], S[J + 1]), max(S[J], S[J + 1]))

a5 = arctan_enclosure(Fraction(1, 5))
a239 = arctan_enclosure(Fraction(1, 239))
PI = (16 * a5[0] - 4 * a239[1], 16 * a5[1] - 4 * a239[0])
assert Fraction("3.14159265358979323846264338327950288419716939937510") < PI[0] < PI[1] < Fraction("3.14159265358979323846264338327950288419716939937511")

def zeta_even(n2):
    """Enclosure of zeta(n2) for even n2 >= 2."""
    c = abs(B[n2]) * 2 ** n2 / (2 * factorial(n2))
    return (c * PI[0] ** n2, c * PI[1] ** n2)

# gamma enclosure from H_D (D = 10^4, m = 6), ln p via atanh series --- same derivation as route (b)
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
DG, MG = 10 ** 4, 6
N, L = 0, 1
spf = list(range(DG + 1))
for i in range(2, int(DG ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, DG + 1, i):
            if spf[j] == j:
                spf[j] = i
for d in range(1, DG + 1):
    if d > 1:
        q, x = spf[d], d
        while x % q == 0:
            x //= q
        if x == 1:
            L *= q; N *= q
    N += L // d
HD = Fraction(N, L)
A_D = 1 / (2 * Fraction(DG)) - sum(B[2 * j] / (2 * j * Fraction(DG) ** (2 * j)) for j in range(1, MG + 1))
RD = abs(B[2 * MG + 2]) / (2 * MG + 2) / Fraction(DG) ** (2 * MG + 2)
ln5 = ln_int_enclosure(5)
lnD = (4 * (LN2[0] + ln5[0]), 4 * (LN2[1] + ln5[1]))
GAMMA = (HD - A_D - lnD[1] - RD, HD - A_D - lnD[0] + RD)
assert Fraction("0.5772156649015328606065120900824024310421") < GAMMA[0] < GAMMA[1] < Fraction("0.5772156649015328606065120900824024310422")

# primes
LIM = 2000
sieve = bytearray([1]) * (LIM + 1); sieve[0] = sieve[1] = 0
for i in range(2, int(LIM ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
primes = [i for i in range(LIM + 1) if sieve[i]]

def fmt(fr, digits=40):
    return "%.*e" % (12, float(fr))

results = {}
print("uniform bound UB(k) >= G(k), UB strictly decreasing in k (exact rationals; pi, gamma enclosed)")
for m in (3, 4, 5, 6, 7):
    coef = {1: Fraction(1, 2)}
    for j in range(1, m + 1):
        coef[2 * j] = -B[2 * j] / (2 * j)
    C_m = abs(B[2 * m + 2]) / (2 * m + 2)
    A1 = sum(coef.values())                                   # A_m(1) = 1/2 - sum B_{2j}/(2j)
    # constant part, upper bound:  sum_{c_s<0} c_s / zeta(s)  ->  need upper bound of zeta(s) (c_s<0)
    const_hi = Fraction(0)
    for s, c in coef.items():
        if c < 0:
            const_hi += c / zeta_even(s)[1]
    const_hi += 1 - GAMMA[0] - A1
    z1 = zeta_even(2 * m + 2); z2 = zeta_even(4 * m + 4)
    tail = C_m * (z1[1] / z2[0] - 1)
    const_hi += tail
    P = {s: Fraction(1) for s in coef if coef[s] > 0}
    first = None
    table = {}
    for k in range(1, 130):
        p = primes[k - 1]
        for s in P:
            P[s] *= 1 - Fraction(1, p ** s)
        if k < 2:
            continue
        UB = sum(coef[s] * P[s] for s in P) + const_hi
        if 86 <= k <= 96:
            table[k] = "%.6e" % float(UB)
        if UB < 0 and first is None:
            first = k
    results[m] = {"first_k_with_UB_negative": first, "p_k0": primes[first - 1], "UB_table_86_96": table,
                  "tail_constant_C_m_times_(zeta(2m+2)/zeta(4m+4)-1)": "%.3e" % float(tail),
                  "eps1_minus_A_m(1)_lower": "%.12e" % float(1 - GAMMA[1] - A1)}
    print("  m=%d: eps(1)-A_m(1) = %+.8e, uniform tail bound = %.3e, UB(91) = %+.4e, UB(92) = %+.4e  ->  UB(k) < 0 for all k >= %d (p_k0 = %d)"
          % (m, float(1 - GAMMA[1] - A1), float(tail), float(Fraction(table[91])), float(Fraction(table[92])), first, primes[first - 1]))

k0 = min(r["first_k_with_UB_negative"] for r in results.values())
print("CONCLUSION: G(k) <= UB(k) < 0 for every k >= %d, hence sign a(p_k#) = (-1)^{k+1} for all k >= %d." % (k0, k0))
out = {"statement": "for all k >= k0: G(k) = (-1)^k F(p_k#) < 0, i.e. a(p_k#) < 0 for even k >= k0 and a(p_k#) > 0 for odd k >= k0",
       "k0": k0,
       "pi_enclosure": ["%.60f" % float(PI[0]), "%.60f" % float(PI[1])],
       "pi_enclosure_width": float(PI[1] - PI[0]),
       "gamma_enclosure_width": float(GAMMA[1] - GAMMA[0]),
       "per_m": results,
       "justification": __doc__}
with open(os.path.join(HERE, "results_uniform_bound.json"), "w") as f:
    json.dump(out, f, indent=1)
print("wrote results_uniform_bound.json")
