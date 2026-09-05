# Minimal rigorous interval arithmetic on top of MPFR (gmpy2) with directed rounding.
# Independent of the proof's tools (fractions + decimal table, mpmath.iv). Every endpoint is produced by an
# MPFR operation rounded toward -inf (lo) or +inf (hi); MPFR guarantees correct rounding for +,-,*,/,log,const_pi,const_euler.
import gmpy2
from gmpy2 import mpfr, mpq
from fractions import Fraction
PREC = 160
CD = gmpy2.context(precision=PREC, round=gmpy2.RoundDown)
CU = gmpy2.context(precision=PREC, round=gmpy2.RoundUp)
def dn(f):
    gmpy2.set_context(CD); return f()
def up(f):
    gmpy2.set_context(CU); return f()
class Iv:
    __slots__ = ('lo', 'hi')
    def __init__(self, lo, hi):
        if not lo <= hi: raise ValueError(("bad interval", lo, hi))
        self.lo, self.hi = lo, hi
    @staticmethod
    def of(x):
        if isinstance(x, Iv): return x
        if isinstance(x, Fraction): x = mpq(x.numerator, x.denominator)
        return Iv(dn(lambda: mpfr(x)), up(lambda: mpfr(x)))
    def __add__(a, b):
        b = Iv.of(b); return Iv(dn(lambda: a.lo + b.lo), up(lambda: a.hi + b.hi))
    __radd__ = __add__
    def __neg__(a): return Iv(-a.hi, -a.lo)
    def __sub__(a, b):
        b = Iv.of(b); return Iv(dn(lambda: a.lo - b.hi), up(lambda: a.hi - b.lo))
    def __rsub__(a, b): return Iv.of(b) - a
    def __mul__(a, b):
        b = Iv.of(b)
        c = [(x, y) for x in (a.lo, a.hi) for y in (b.lo, b.hi)]
        return Iv(min(dn(lambda: x*y) for x, y in c), max(up(lambda: x*y) for x, y in c))
    __rmul__ = __mul__
    def recip(a):
        if not a.lo > 0: raise ValueError("recip of interval not strictly positive")
        return Iv(dn(lambda: 1/a.hi), up(lambda: 1/a.lo))
    def __truediv__(a, b): return a * Iv.of(b).recip()
    def __rtruediv__(a, b): return Iv.of(b) * a.recip()
    def __pow__(a, n):
        r = Iv.of(1)
        for _ in range(n): r = r * a
        return r
    def width(a): return up(lambda: a.hi - a.lo)
    def __repr__(a): return "[%.17e, %.17e]" % (float(a.lo), float(a.hi))
    def strfull(a): return "[%s, %s]" % (str(a.lo)[:40], str(a.hi)[:40])
PI = Iv(dn(gmpy2.const_pi), up(gmpy2.const_pi))
GAMMA = Iv(dn(gmpy2.const_euler), up(gmpy2.const_euler))
def LOG(x):  # x exact int
    return Iv(dn(lambda: gmpy2.log(mpfr(x))), up(lambda: gmpy2.log(mpfr(x))))
LN2 = LOG(2)
assert PI.lo < PI.hi and GAMMA.lo < GAMMA.hi and LN2.lo < LN2.hi, "directed rounding of constants not effective"
