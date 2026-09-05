#!/usr/bin/env python3
# Exact-rational checks: (i) a(n) from the defining relation, integrality, OEIS terms; (ii) (B.2) vs exact F(n) for
# squarefree n<=N with omega>=2, and (6.1) for non-squarefree n; (iii) brute-force G(p_k#) over all 2^k divisors, k<=12.
from fractions import Fraction as Fr
from math import factorial
from sympy import factorint, divisors, mobius, primerange
from mpmath import mp, mpf, euler, harmonic, fprod, zeta
mp.dps = 50
N = 4000
H = [Fr(0)]
for i in range(1, N+1): H.append(H[-1] + Fr(1, i))
# (i) a(n) from sum_{d|n} a(d)/d! = H_n, n<=40
a = {}
for n in range(1, 41):
    a[n] = factorial(n) * (H[n] - sum(a[d]/factorial(d) for d in divisors(n) if d < n))
print("a(1..12) =", [a[n] for n in range(1, 13)]); print("all integers n<=40:", all(a[n].denominator == 1 for n in a))
print("a(30) =", a[30], " (OEIS: -22690644647302814715858124800000)")
# F(n) via Moebius inversion equals a(n)/n! :
def F(n): return sum(int(mobius(n//d)) * H[d] for d in divisors(n))
print("Moebius inversion consistent n<=40:", all(F(n) == a[n]/factorial(n) for n in a))
omega = lambda n: len(factorint(n))
viol = [n for n in range(1, 41) if (a[n] < 0) != (omega(n) % 2 == 1 and omega(n) >= 3)]
print("violations of conjecture for n<=40:", viol)
# (ii)
c4 = mpf(2897)/5040 - euler
def P(ps, s): return fprod([1 - mpf(p)**(-s) for p in ps])
worst_sq = 0; cnt_sq = 0; worst_ns = 0; cnt_ns = 0; bad = []
for n in range(2, N+1):
    f = factorint(n); ps = sorted(f); k = len(ps)
    if k < 2: continue
    q = 1
    for p, e in f.items(): q *= p**(e-1)
    Gx = (-1)**k * F(n); Gx = mpf(Gx.numerator)/Gx.denominator
    if q == 1:
        Gm = P(ps,1)/2 - P(ps,2)/12 + P(ps,4)/120 - P(ps,6)/252 + P(ps,8)/240 + c4
        Eb = (fprod([1 + mpf(p)**(-10) for p in ps]) - 1)/132
        cnt_sq += 1; r = abs(Gx-Gm)/Eb; worst_sq = max(worst_sq, r)
        if r > 1: bad.append(n)
        assert Eb < mpf(1)/132000
    else:
        Gm = P(ps,1)/(2*q) - P(ps,2)/(12*q**2) + P(ps,4)/(120*q**4) - P(ps,6)/(252*q**6) + P(ps,8)/(240*q**8)
        Eb = fprod([1 + mpf(p)**(-10) for p in ps])/(132*mpf(q)**10)
        cnt_ns += 1; r = abs(Gx-Gm)/Eb; worst_ns = max(worst_ns, r)
        if r > 1: bad.append(n)
    if Gx <= 0: bad.append(('G<=0', n))
print(f"(B.2): {cnt_sq} squarefree n<={N}, omega>=2: max |G_exact-G_main|/E_bound = {float(worst_sq):.5f}")
print(f"(6.1): {cnt_ns} non-squarefree n<={N}, omega>=2: max ratio = {float(worst_ns):.5f}")
print("violations of bound or G<=0:", bad)
# (iii) brute force primorials
primes = list(primerange(2, 40))
for k in range(2, 13):
    ps = primes[:k]; divs = [(1, 0)]
    for p in ps: divs += [(d*p, c+1) for d, c in divs]
    Gb = (-1)**k * sum((-1)**(k-c) * harmonic(d) for d, c in divs)
    Gm = P(ps,1)/2 - P(ps,2)/12 + P(ps,4)/120 - P(ps,6)/252 + P(ps,8)/240 + c4
    Eb = (fprod([1 + mpf(p)**(-10) for p in ps]) - 1)/132
    print(f"k={k:2d} brute G={float(Gb):+.15e} closed={float(Gm):+.15e} |diff|={float(abs(Gb-Gm)):.2e} <= E_bound={float(Eb):.2e}: {abs(Gb-Gm) <= Eb}")
print("zeta(10)-1 =", zeta(10)-1)
