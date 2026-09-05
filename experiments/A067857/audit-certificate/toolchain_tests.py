#!/usr/bin/env python3
"""Sanity tests of the tools the certificate relies on, and of my own planned arithmetic.
(1) mpmath.iv: basic ops round outward?  (2) gmpy2 mpfr: rounding modes honoured?  (3) certificate's decimal
outward-rounding helpers (copied verbatim from route_a_iv.py / route_b_exact.py) on negative inputs."""
from fractions import Fraction
from mpmath import iv, mp
from mpmath.libmp import to_rational
import gmpy2
from gmpy2 import mpfr, mpq

iv.prec = 64
def ends(x):
    a, b = x._mpi_
    return Fraction(*to_rational(a)), Fraction(*to_rational(b))
# (1) mpmath.iv
x = iv.mpf(1) / iv.mpf(3)
lo, hi = ends(x)
print("iv 1/3 :", lo < Fraction(1, 3) < hi, "width", float(hi - lo))
y = iv.mpf(7) ** 10
lo, hi = ends(y); print("iv 7**10 exact int:", lo == 7 ** 10 == hi)
z = iv.mpf(1) / (iv.mpf(479) ** 9)
lo, hi = ends(z); print("iv 1/479^9 :", lo < Fraction(1, 479 ** 9) < hi)
w = iv.mpf(1) / iv.mpf(3) - iv.mpf(1) / iv.mpf(7)
lo, hi = ends(w); print("iv 1/3-1/7 :", lo < Fraction(1, 3) - Fraction(1, 7) < hi)
h = iv.mpf([0, iv.mpf(1) / iv.mpf(3)])
lo, hi = ends(h); print("iv hull [0, 1/3]:", lo == 0 and hi > Fraction(1, 3) - Fraction(1, 10**18))
# accumulate 10^4 terms of 1/i and compare with exact H
H = iv.mpf(0); He = Fraction(0)
for i in range(1, 10001):
    H = H + iv.mpf(1) / iv.mpf(i); He += Fraction(1, i)
lo, hi = ends(H); print("iv H_10000 encloses exact:", lo <= He <= hi, "width", float(hi - lo))
# (2) gmpy2 directed rounding
ctx = gmpy2.get_context(); ctx.precision = 100
ctx.round = gmpy2.RoundDown; d = mpfr(1) / 3
ctx.round = gmpy2.RoundUp;   u = mpfr(1) / 3
print("gmpy2 1/3 down<1/3<up:", mpq(d) < Fraction(1, 3) < mpq(u) if False else (Fraction(int(mpq(d).numerator), int(mpq(d).denominator)) < Fraction(1, 3) < Fraction(int(mpq(u).numerator), int(mpq(u).denominator))))
print("gmpy2 exact ratio:", Fraction(*d.as_integer_ratio()) < Fraction(1,3) < Fraction(*u.as_integer_ratio()))
ctx.round = gmpy2.RoundDown; ld = gmpy2.log(mpfr(2))
ctx.round = gmpy2.RoundUp;   lu = gmpy2.log(mpfr(2))
print("gmpy2 log2 down/up differ:", ld < lu, float(lu - ld))
# (3) certificate helper functions (verbatim copies)
def fmt_fixed(q, digits):
    sign = '-' if q < 0 else ''
    s = str(abs(q)).rjust(digits + 1, '0')
    return sign + s[:-digits] + '.' + s[-digits:]
def dec_floor(fr, digits):
    return fmt_fixed((fr.numerator * 10 ** digits) // fr.denominator, digits)
def dec_ceil(fr, digits):
    return fmt_fixed(-((-fr.numerator * 10 ** digits) // fr.denominator), digits)
import random
random.seed(1)
bad = 0
for _ in range(20000):
    n = random.randint(-10**12, 10**12); d = random.randint(1, 10**9); dg = random.randint(1, 30)
    fr = Fraction(n, d)
    lo = Fraction(dec_floor(fr, dg)); hi = Fraction(dec_ceil(fr, dg))
    if not (lo <= fr <= hi and hi - lo <= Fraction(1, 10 ** dg)):
        bad += 1
print("certificate dec_floor/dec_ceil outward on 20000 random fractions (incl. negative): failures =", bad)
print(dec_floor(Fraction(-1, 3), 5), dec_ceil(Fraction(-1, 3), 5), dec_floor(Fraction(-5, 1), 3), dec_ceil(Fraction(-5, 1), 3), dec_floor(Fraction(0), 3))
