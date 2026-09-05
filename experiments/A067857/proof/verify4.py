from fractions import Fraction
from sympy import primerange
from mpmath import mp, mpf, iv, euler, pi, zeta, harmonic, log, bernoulli
mp.dps = 40
primes = list(primerange(2, 150000))
assert primes[91] == 479
# (a) interval arithmetic UB(92)
iv.dps = 30
def ivprod(s, k):
    r = iv.mpf(1)
    for p in primes[:k]: r *= (1 - iv.mpf(1)/iv.mpf(p)**s)
    return r
UBiv = ivprod(1,92)/2 + ivprod(4,92)/120 + ivprod(8,92)/240 - 1/(2*iv.pi**2) - iv.mpf(15)/(4*iv.pi**6) + iv.mpf(2897)/5040 - iv.euler + iv.mpf(1)/132000
print("(a) interval UB(92) =", UBiv)
# (b) coarse certificate: pi<=3.1415926536, gamma>=0.5772156649, products rounded UP at 10 decimals
def Pex(s, k):
    r = Fraction(1)
    for p in primes[:k]: r *= 1 - Fraction(1, p**s)
    return r
def round_up(fr, dec=10):
    q = 10**dec
    return Fraction(-((-fr.numerator*q)//fr.denominator), q)   # ceil
def round_down(fr, dec=10):
    q = 10**dec
    return Fraction((fr.numerator*q)//fr.denominator, q)        # floor
P1, P4, P8 = Pex(1,92), Pex(4,92), Pex(8,92)
print("    exact P1(92) = %s / %s (approx %s)" % (str(P1.numerator)[:20]+'...', str(P1.denominator)[:20]+'...', mpf(P1.numerator)/P1.denominator))
pi_hi = Fraction('3.1415926536'); g_lo = Fraction('0.5772156649')
terms = [("(1/2)P1(92)", round_up(P1/2)), ("(1/120)P4(92)", round_up(P4/120)), ("(1/240)P8(92)", round_up(P8/240)),
         ("-1/(2 pi^2)", round_up(-1/(2*pi_hi**2))), ("-15/(4 pi^6)", round_up(-Fraction(15,4)/pi_hi**6)),
         ("2897/5040", round_up(Fraction(2897,5040))), ("-gamma", round_up(-g_lo)), ("1/132000", round_up(Fraction(1,132000)))]
tot = Fraction(0)
for name, v in terms:
    print(f"    {name:14s} <= {float(v):+.10f}"); tot += v
print("    sum of upper bounds =", float(tot), " (<0:", tot < 0, ")")
# same with 12-decimal reference values (not rounded) for the table
mp.dps = 30
ref = [("(1/2)P1(92)", mpf(P1.numerator)/P1.denominator/2), ("(1/120)P4(92)", mpf(P4.numerator)/P4.denominator/120), ("(1/240)P8(92)", mpf(P8.numerator)/P8.denominator/240),
       ("-1/(2 pi^2)", -1/(2*pi**2)), ("-15/(4 pi^6)", -mpf(15)/(4*pi**6)), ("2897/5040", mpf(2897)/5040), ("-gamma", -euler), ("1/132000", mpf(1)/132000)]
print("    reference values (12 decimals):"); s = mpf(0)
for name, v in ref: print(f"    {name:14s} = {v}"); s += v
print("    UB(92) =", s)
# UB(91) reference
s91 = mpf(Pex(1,91).numerator)/Pex(1,91).denominator/2 + mpf(Pex(4,91).numerator)/Pex(4,91).denominator/120 + mpf(Pex(8,91).numerator)/Pex(8,91).denominator/240 - 1/(2*pi**2) - mpf(15)/(4*pi**6) + mpf(2897)/5040 - euler + mpf(1)/132000
print("    UB(91) =", s91)
# (c) k6, k12
def Q(k):
    r = mpf(1)
    for p in primes[:k]: r *= 1 + mpf(1)/p
    return r
k6 = next(k for k in range(1, 2000) if Q(k) > 6); k12 = next(k for k in range(1, 20000) if Q(k) > 12)
print(f"(c) least k with Q(p_k#)>6: k6={k6}, p={primes[k6-1]}, Q(k6-1)={Q(k6-1)}, Q(k6)={Q(k6)}")
print(f"    least k with Q(p_k#)>12: k12={k12}, p={primes[k12-1]}, Q(k12-1)={Q(k12-1)}, Q(k12)={Q(k12)}")
# (e) positivity constants (Prop 7)
c4 = mpf(2897)/5040 - euler
br_sq = 1/(120*zeta(4)) - mpf(1)/252 + 1/(240*zeta(8)) + c4 - mpf(1)/132000
print("(e) squarefree bracket lower bound 1/(120 z4) - 1/252 + 1/(240 z8) + c4 - 1/132000 =", br_sq)
q = 2
br_q = 1/(120*zeta(4)) - mpf(1)/(252*q**2) - mpf('1.001')/(132*q**6)
print("    q>=2 bracket lower bound (times q^-4): 1/(120 z4) - 1/(252*4) - 1.001/(132*64) =", br_q)
# (g) heuristic bracket range over k
def Pk(s,k):
    r = mpf(1)
    for p in primes[:k]: r *= 1 - mpf(p)**(-s)
    return r
vals = [Pk(4,k)/120 - Pk(6,k)/252 + Pk(8,k)/240 + c4 for k in range(1, 400)]
print("(g) bracket (1/120)P4-(1/252)P6+(1/240)P8+c4 over k=1..399: min", min(vals), "max", max(vals), " limit", 1/(120*zeta(4)) - 1/(252*zeta(6)) + 1/(240*zeta(8)) + c4)
# (d) non-squarefree n_k = 2*p_k#, q=2, with R_4(2) exact
R42 = harmonic(2) - log(2) - euler - (mpf(1)/4 - sum(bernoulli(2*j)/(2*j*mpf(2)**(2*j)) for j in range(1,5)))
print("(d) R_4(2) = eps(2) - A_4(2) =", R42)
Ecap2 = mpf(1)/(132*1024) * mpf('1e-3')    # (1/132) 2^-10 (prod(1+p^-10)-1) < (1/132) 2^-10 (zeta(10)-1) < ... 
print("    residual error cap (1/132)2^-10*1e-3 =", Ecap2)
def G2(k):
    return Pk(1,k)/4 - Pk(2,k)/48 + Pk(4,k)/1920 - Pk(6,k)/16128 + Pk(8,k)/61440 + R42
lo, hi = 9000, 9400
while hi - lo > 1:
    mid = (lo+hi)//2
    if G2(mid) < 0: hi = mid
    else: lo = mid
print(f"    G(2*p_k#): k={lo}: {G2(lo)} ; k={hi}: {G2(hi)} ; p_{hi} = {primes[hi-1]}; certified: {G2(lo) > Ecap2 and G2(hi) < -Ecap2}")
def UB2(k):
    return Pk(1,k)/4 + Pk(4,k)/1920 + Pk(8,k)/61440 - (6/pi**2)/48 - (945/pi**6)/16128 + R42 + Ecap2
k1 = next(k for k in range(hi, hi+2000) if UB2(k) < 0)
print(f"    least k with UB2(k)<0: {k1}, p_k={primes[k1-1]}, UB2(k1-1)={UB2(k1-1)}, UB2(k1)={UB2(k1)}")
# positivity of G(2 p_k#) for all k < hi (lower bound)
print("    min over k<hi of G2(k)-Ecap2 > 0 ?", min(G2(k) - Ecap2 for k in range(2, hi)) > 0)
# digits of 2*p_hi#
import math
print("    log10(2*p_hi#) ~", sum(math.log10(p) for p in primes[:hi]) + math.log10(2))
