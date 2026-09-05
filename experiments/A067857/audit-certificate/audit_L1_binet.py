#!/usr/bin/env python3
"""Audit of lemma L1 (the only non-elementary input of the certificate).

L1: for real x>0, integer m>=0:  psi(x) = ln x - 1/(2x) - sum_{j=1}^m B_{2j}/(2j x^{2j}) + R_m(x),
    sign R_m(x) = sign(-B_{2m+2}),  |R_m(x)| < |B_{2m+2}| / ((2m+2) x^{2m+2}).
My own proof (written out in README.md of this directory) starts from Binet's second formula for psi,
    psi(x) = ln x - 1/(2x) - 2 int_0^oo t / ((t^2+x^2)(e^{2 pi t}-1)) dt            (x>0)
and the moments  int_0^oo t^{2j+1}/(e^{2 pi t}-1) dt = |B_{2j+2}| / (4(j+1)).
This script checks numerically, at 60-100 digits: (a) Binet's formula, (b) the moments, (c) L1 sign and magnitude
for many (x, m) including x = 1 (the smallest d used) and the D-values used by the certificates and by me,
(d) H_d = psi(d) + 1/d + gamma, (e) the identity Ein(x) = gamma + ln x + E_1(x) that I use for gamma.
These are numerical sanity checks; the proof is the derivation in README.md."""
from mpmath import mp, mpf, psi, log, quad, inf, exp, pi, bernoulli, harmonic, euler, e1, factorial
mp.dps = 60
print("(a) Binet's second formula for psi, |lhs - rhs| at 60 digits:")
for x in [mpf('0.5'), mpf(1), mpf(2), mpf(3), mpf(10), mpf(100), mpf(3000)]:
    f = lambda t: t / ((t * t + x * x) * (exp(2 * pi * t) - 1))
    integral = quad(f, [0, 1, 10, inf])
    rhs = log(x) - 1 / (2 * x) - 2 * integral
    print("   x=%-6s |psi(x) - rhs| = %.2e" % (mp.nstr(x, 5), abs(psi(0, x) - rhs)))
print("(b) moments int_0^oo t^{2j+1}/(e^{2 pi t}-1) dt vs |B_{2j+2}|/(4(j+1)):")
for j in range(0, 8):
    f = lambda t: t ** (2 * j + 1) / (exp(2 * pi * t) - 1)
    integral = quad(f, [0, 1, 5, 20, inf])
    target = abs(bernoulli(2 * j + 2)) / (4 * (j + 1))
    print("   j=%d  |diff| = %.2e" % (j, abs(integral - target)))
mp.dps = 110
print("(c) L1 sign and magnitude, 110 digits (ratio = |R_m(x)| / bound must be < 1, sign must match -B_{2m+2}):")
worst = 0; bad = 0; cnt = 0
xs = [1, 2, 3, 4, 5, 6, 7, 10, 30, 100, 1000, 3000, 3001, 10**4, 10**5, 65536, mpf('0.5'), mpf('1.5'), mpf('0.1')]
for x in xs:
    x = mpf(x)
    for m in range(0, 13):
        S = log(x) - 1 / (2 * x) - sum(bernoulli(2 * j) / (2 * j * x ** (2 * j)) for j in range(1, m + 1))
        R = psi(0, x) - S
        bound = abs(bernoulli(2 * m + 2)) / ((2 * m + 2) * x ** (2 * m + 2))
        sgn_expected = -1 if bernoulli(2 * m + 2) > 0 else 1
        ratio = abs(R) / bound
        cnt += 1
        ok = (R * sgn_expected > 0) and ratio < 1
        if not ok:
            bad += 1; print("   FAIL x=%s m=%d R=%s bound=%s" % (x, m, R, bound))
        worst = max(worst, ratio)
print("   %d/%d (x,m) cases pass; largest ratio |R|/bound = %.6f" % (cnt - bad, cnt, worst))
print("(d) H_d - (psi(d) + 1/d + gamma):", max(abs(harmonic(d) - (psi(0, mpf(d)) + mpf(1) / d + euler)) for d in [1, 2, 3, 7, 3000, 65536]))
print("(e) Ein(x) - ln x - E_1(x) - gamma, Ein by its power series (A&S 5.1.11 / DLMF 6.2.3-6.2.4):")
for x in [1, 10, 256]:
    x = mpf(x)
    Ein = mp.nsum(lambda k: (-1) ** (k - 1) * x ** k / (k * factorial(k)), [1, inf])
    print("   x=%s  |Ein - ln x - E_1(x) - gamma| = %.2e" % (mp.nstr(x, 4), abs(Ein - log(x) - e1(x) - euler)))
print("(f) e > 2.7 (used for E_1(256) < e^{-256}/256 < 2.7^{-256}/256):", exp(1) > mpf('2.7'), "  2.7^-256/256 = %.3e" % (mpf('2.7') ** -256 / 256))
