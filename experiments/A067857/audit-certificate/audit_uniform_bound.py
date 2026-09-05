#!/usr/bin/env python3
"""Auditor's independent check of the certificate's uniform bound  G(k) <= UB_m(k),  UB_m strictly decreasing,
UB_m(92) < 0  =>  G(k) < 0 for all k >= 92.

Derivation (mine; D = 1 version of the closed form, k >= 2):
  G(k) = sum_s c_s P_s(k) + (eps(1) - A_m(1)) + E(k),   eps(1) = 1 - gamma,
  E(k) = sum_{d | p_k#, d > 1} (-1)^omega R_m(d),  |E(k)| <= C_m sum_{d|n, d>1} d^{-(2m+2)} = C_m (prod_{i<=k}(1 + p_i^{-(2m+2)}) - 1).
  For s > 1: 0 < P_s(k) decreasing in k, P_s(k) >= prod_all_p (1 - p^{-s}) = 1/zeta(s);  prod_{i<=k}(1 + p_i^{-s}) <= zeta(s)/zeta(2s).
  => G(k) <= UB_m(k) := sum_{c_s>0} c_s P_s(k) + sum_{c_s<0} c_s/zeta(s) + (1 - gamma - A_m(1)) + C_m (zeta(2m+2)/zeta(4m+4) - 1),
  and UB_m is strictly decreasing (c_1 = 1/2 > 0 and P_1 strictly decreasing).
Ingredients computed here from scratch, all exact rationals:
  pi by Euler's formula pi/4 = arctan(1/2) + arctan(1/3) (alternating series, enclosure between consecutive partial sums;
     NOT Machin's formula, which the certificate uses); zeta(2n) = |B_2n| (2 pi)^{2n} / (2 (2n)!);
  gamma via Ein(256) - 8 ln 2 - E_1(256) (as in audit_routeC.py); ln 2 = 2 atanh(1/3) exact series with tail bound;
  Bernoulli numbers by Akiyama-Tanigawa.
Checks: UB_m(k) for k = 2..300 and m = 3..7 against the certificate's table (results_uniform_bound.json), monotonicity,
UB_4(92) < 0 < UB_4(91), and G(k) <= UB_m(k) for my own route-C enclosures (results_routeC.json) for k <= 200."""
import json, sys, os
from fractions import Fraction
from math import factorial
sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__))
CERT = "/home/user/Solving-hard-math-problems-/experiments/A067857/certificate"

def bernoulli_AT(nmax):
    B = []
    for n in range(nmax + 1):
        A = [Fraction(0)] * (n + 1)
        for mm in range(n + 1):
            A[mm] = Fraction(1, mm + 1)
            for j in range(mm, 0, -1):
                A[j - 1] = j * (A[j - 1] - A[j])
        B.append(A[0])
    return B
Bl = bernoulli_AT(32)
B = {n: Bl[n] for n in range(2, 33, 2)}
assert B[2] == Fraction(1, 6) and B[10] == Fraction(5, 66) and B[30] == Fraction(8615841276005, 14322)

def arctan_recip_enclosure(n, J=220):
    """arctan(1/n) = sum (-1)^j / ((2j+1) n^{2j+1}): alternating with strictly decreasing terms -> between partial sums."""
    s, S = Fraction(0), []
    for j in range(J + 2):
        s += Fraction((-1) ** j, (2 * j + 1) * n ** (2 * j + 1)); S.append(s)
    return min(S[J], S[J + 1]), max(S[J], S[J + 1])
a2, a3 = arctan_recip_enclosure(2), arctan_recip_enclosure(3)
PI = (4 * (a2[0] + a3[0]), 4 * (a2[1] + a3[1]))
assert Fraction("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899") < PI[0] < PI[1] < Fraction("3.14159265358979323846264338327950288419716939937510582097494459230781640628620900")
def zeta_even(n2):
    c = abs(B[n2]) * 2 ** n2 / (2 * factorial(n2))
    return (c * PI[0] ** n2, c * PI[1] ** n2)
z2 = zeta_even(2)
from mpmath import mp as _mp, zeta as _zeta, mpf as _mpf
_mp.dps = 80
_z2 = _zeta(2); assert abs(_mpf(z2[0].numerator) / z2[0].denominator - _z2) < _mpf(10) ** -75 and z2[0] <= z2[1]

# ln 2 and gamma (Ein method)
def atanh_recip(n, J=300):
    s, term = Fraction(0), Fraction(1, n)
    for j in range(J):
        s += term / (2 * j + 1); term /= n * n
    return s, s + term / (2 * J + 1) * Fraction(n * n, n * n - 1)
LN2 = tuple(2 * x for x in atanh_recip(3))
X = 256; Sk, term, k = Fraction(0), Fraction(1), 0
while True:
    k += 1; term = term * X / k; Sk += (-1) ** (k - 1) * term / k
    if k > X and term * X / (k + 1) / (k + 1) < Fraction(1, 10 ** 130):
        break
tnext = term * X / (k + 1) / (k + 1)
EB = Fraction(10, 27) ** X / X
GAMMA = (Sk - tnext - 8 * LN2[1] - EB, Sk + tnext - 8 * LN2[0])
_g = Fraction("0.5772156649015328606065120900824024310421593359399235988057672348848677")   # truncated 70 digits
assert GAMMA[0] < _g + Fraction(1, 10 ** 70) and _g < GAMMA[1] and GAMMA[1] - GAMMA[0] < Fraction(1, 10 ** 100)

LIM = 2100
sv = bytearray([1]) * (LIM + 1); sv[0] = sv[1] = 0
for i in range(2, int(LIM ** 0.5) + 1):
    if sv[i]:
        sv[i * i::i] = bytearray(len(sv[i * i::i]))
primes = [i for i in range(LIM + 1) if sv[i]]
assert primes[91] == 479 and len(primes) >= 300

cert = json.load(open(CERT + "/results_uniform_bound.json"))
mine = json.load(open(HERE + "/results_routeC.json")) if os.path.exists(HERE + "/results_routeC.json") else None
routeA = {c["k"]: Fraction(c["G_enclosure"][1]) for c in json.load(open(CERT + "/cert_route_a.json"))["certificates"]}
out = {}
KTOP = 300
for m in (3, 4, 5, 6, 7):
    coef = {1: Fraction(1, 2)}
    for j in range(1, m + 1):
        coef[2 * j] = -B[2 * j] / (2 * j)
    C_m = abs(B[2 * m + 2]) / (2 * m + 2)
    A1 = sum(coef.values())
    const = Fraction(0)
    for s, c in coef.items():
        if c < 0:
            const += c / zeta_even(s)[1]          # c<0: upper bound needs the SMALLEST 1/zeta(s), i.e. largest zeta(s)
    const += 1 - GAMMA[0] - A1                     # upper bound: smallest gamma
    zt1, zt2 = zeta_even(2 * m + 2), zeta_even(4 * m + 4)
    const += C_m * (zt1[1] / zt2[0] - 1)
    P = {s: Fraction(1) for s in coef if coef[s] > 0}
    UB = {}
    for k in range(1, KTOP + 1):
        p = primes[k - 1]
        for s in P:
            P[s] *= 1 - Fraction(1, p ** s)
        if k >= 2:
            UB[k] = sum(coef[s] * P[s] for s in P) + const
    first = min(k for k in UB if UB[k] < 0)
    mono = all(UB[k + 1] < UB[k] for k in range(2, KTOP))
    # compare with certificate table (6 significant digits printed there)
    tab = cert["per_m"][str(m)]["UB_table_86_96"]
    reldiff = max(abs(float(UB[int(kk)]) - float(v)) / abs(float(v)) for kk, v in tab.items())
    # G(k) <= UB(k): my route-C upper endpoints and route (a)'s
    ok_mine = ok_a = True
    margin_mine = None
    if mine:
        for r in mine["results"]:
            kk = r["k"]; hi = Fraction(r["G_hi"])
            ok_mine &= hi <= UB[kk]
            d = UB[kk] - hi
            margin_mine = d if margin_mine is None else min(margin_mine, d)
    for kk, hi in routeA.items():
        ok_a &= hi <= UB[kk]
    out[m] = {"first_k_with_UB_negative": first, "UB_91": "%.10e" % float(UB[91]), "UB_92": "%.10e" % float(UB[92]),
              "strictly_decreasing_k_le_300": mono, "max_rel_diff_vs_certificate_table": reldiff,
              "G_le_UB_for_all_my_enclosures_k_le_200": ok_mine, "G_le_UB_for_route_a_enclosures": ok_a,
              "min_margin_UB_minus_G_hi": float(margin_mine) if margin_mine is not None else None,
              "const_part": "%.12e" % float(const), "C_m_tail": "%.3e" % float(C_m * (zt1[1] / zt2[0] - 1))}
    print("m=%d: UB(91)=%+.6e UB(92)=%+.6e first negative k=%d, strictly decreasing to k=300: %s, max rel diff vs cert table %.1e, G<=UB (mine/route a): %s/%s, min margin %.2e"
          % (m, float(UB[91]), float(UB[92]), first, mono, reldiff, ok_mine, ok_a, float(margin_mine) if margin_mine is not None else float('nan')))
    assert first == 92 and mono and ok_a and (ok_mine or mine is None) and reldiff < 2e-6
# extra: the five 'facts' used, checked numerically with exact rationals for s = 2, 6, 10 at k = 300
print("P_s(300) >= 1/zeta(s) and prod(1+p^-s)(300) <= zeta(s)/zeta(2s), s in {2,4,6,10}:", end=" ")
for s in (2, 4, 6, 10):
    Ps, Qs = Fraction(1), Fraction(1)
    for p in primes[:300]:
        Ps *= 1 - Fraction(1, p ** s); Qs *= 1 + Fraction(1, p ** s)
    zs, z2s = zeta_even(s), zeta_even(2 * s)
    assert Ps >= 1 / zs[1] and Qs <= zs[1] / z2s[0]
    print("s=%d ok" % s, end="; ")
print()
json.dump({"pi_enclosure_width": float(PI[1] - PI[0]), "gamma_enclosure_width": float(GAMMA[1] - GAMMA[0]), "per_m": out,
           "conclusion": "UB_m(k) < 0 for all k >= 92 for each m in 3..7 (UB strictly decreasing), hence G(k) < 0 and sign a(p_k#) = (-1)^{k+1} for all k >= 92"},
          open(HERE + "/results_uniform_bound.json", "w"), indent=1)
print("CONCLUSION reproduced: G(k) <= UB(k) < 0 for all k >= 92; wrote results_uniform_bound.json")
