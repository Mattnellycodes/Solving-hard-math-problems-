# Part 4: certificate arithmetic for UB(92) < 0 with exact rationals and directed rounding; also interval arithmetic (mpmath.iv).
from fractions import Fraction
from sympy import primerange
from mpmath import mp, mpf, iv, euler, pi, fprod
mp.dps = 50
primes = list(primerange(2, 480)); assert len(primes) == 92 and primes[-1] == 479
def Pexact(s): 
    r = Fraction(1)
    for p in primes: r *= 1 - Fraction(1, p**s)
    return r
P1, P4, P8 = Pexact(1), Pexact(4), Pexact(8)
def dec(fr, digits=15):  # decimal string (truncated)
    return mpf(fr.numerator)/fr.denominator
print("P1(92) =", dec(P1)); print("P4(92) =", dec(P4)); print("P8(92) =", dec(P8))
# rational enclosures: pi in (3.14159265358979, 3.14159265358980), gamma in (0.57721566490153, 0.57721566490154)
pi_lo, pi_hi = Fraction('3.14159265358979'), Fraction('3.14159265358980')
g_lo, g_hi = Fraction('0.57721566490153'), Fraction('0.57721566490154')
assert pi_lo < Fraction(str(pi)) < pi_hi and g_lo < Fraction(str(euler)) < g_hi
# UB(92) = P1/2 + P4/120 + P8/240 - 1/(2 pi^2) - 15/(4 pi^6) + 2897/5040 - gamma + 1/132000
# upper bound: use pi_hi (makes the negative terms least negative) and g_lo
UB_hi = P1/2 + P4/120 + P8/240 - 1/(2*pi_hi**2) - Fraction(15,4)/pi_hi**6 + Fraction(2897,5040) - g_lo + Fraction(1,132000)
print("rigorous upper bound for UB(92):", float(UB_hi), " <0:", UB_hi < 0)
# breakdown of terms
terms = {"P1/2": P1/2, "P4/120": P4/120, "P8/240": P8/240, "-1/(2pi^2)": -1/(2*pi_hi**2), "-15/(4pi^6)": -Fraction(15,4)/pi_hi**6,
         "2897/5040": Fraction(2897,5040), "-gamma": -g_lo, "1/132000": Fraction(1,132000)}
for k_, v in terms.items(): print(f"  {k_:12s} = {float(v):+.12f}")
# same with mpmath interval arithmetic
iv.dps = 30
Pi1 = fprod([1 - iv.mpf(1)/p for p in primes]); Pi4 = fprod([1 - iv.mpf(1)/p**4 for p in primes]); Pi8 = fprod([1 - iv.mpf(1)/p**8 for p in primes])
UBiv = Pi1/2 + Pi4/120 + Pi8/240 - 1/(2*iv.pi**2) - iv.mpf(15)/(4*iv.pi**6) + iv.mpf(2897)/5040 - iv.euler + iv.mpf(1)/132000
print("interval UB(92) =", UBiv)
# sensitivity: if every input is known to +-1e-8, total error <= ?
print("coefficient sum on inputs (1/2+1/120+1/240 for products; d/dpi terms; 1 for gamma):", float(Fraction(1,2)+Fraction(1,120)+Fraction(1,240)), 
      " |d/dpi[-1/(2pi^2)-15/(4pi^6)]| =", float(1/pi**3 + 15*6/(4*pi**7)))
# also: P_1(91) etc for UB(91) (not needed, but shows UB(91)>0)
primes91 = primes[:91]
def Pe(s, ps):
    r = Fraction(1)
    for p in ps: r *= 1 - Fraction(1, p**s)
    return r
UB91 = Pe(1,primes91)/2 + Pe(4,primes91)/120 + Pe(8,primes91)/240 - 1/(2*pi_lo**2) - Fraction(15,4)/pi_lo**6 + Fraction(2897,5040) - g_hi + Fraction(1,132000)
print("rigorous lower bound for UB(91):", float(UB91), " >0:", UB91 > 0)
# Q(91), Q(92)
Q = lambda ps: fprod([1 + mpf(1)/p for p in ps])
print("Q(91) =", Q(primes91), " Q(92) =", Q(primes), " P1(92)*Q(92) = P2(92) =", dec(Pexact(2)))
