# Independent verification, part 2: uniform bound UB(k), certified enclosures of G(p_k#), k<=91 positivity, non-squarefree.
from mpmath import mp, mpf, euler, pi, fprod, fsum, zeta, bernoulli, harmonic, log
from sympy import primerange
mp.dps = 40
primes = list(primerange(2, 200000))
print("p_91, p_92, p_93 =", primes[90], primes[91], primes[92])
def P(s, k): return fprod([1 - mpf(p)**(-s) for p in primes[:k]])
def Qp(s, k): return fprod([1 + mpf(p)**(-s) for p in primes[:k]])
c4 = mpf(2897)/5040 - euler
print("c4 = 2897/5040 - gamma =", c4, "  check:", mpf(1)/2+mpf(1)/12-mpf(1)/120+mpf(1)/252-mpf(1)/240 - euler)
# elementary bound zeta(10)-1 < 2^-10 + 3^-10 + 3^-9/9
zb = mpf(2)**-10 + mpf(3)**-10 + mpf(3)**-9/9
print("zeta(10)-1 =", zeta(10)-1, " elementary bound:", zb, " < 1e-3:", zb < mpf('1e-3'))
Ecap = mpf(1)/132000
def G_main(k): return P(1,k)/2 - P(2,k)/12 + P(4,k)/120 - P(6,k)/252 + P(8,k)/240 + c4
def E_bound(k): return (Qp(10,k) - 1)/132
def UB(k): return P(1,k)/2 + P(4,k)/120 + P(8,k)/240 - 1/(2*pi**2) - mpf(15)/(4*pi**6) + c4 + Ecap
print("\nk, G_main, E_bound, [G_lo, G_hi], UB(k), Q(k)=prod(1+1/p)")
for k in list(range(86, 100)):
    g = G_main(k); e = E_bound(k)
    print(f"k={k:3d} p_k={primes[k-1]:4d}  G in [{float(g-e):+.6e}, {float(g+e):+.6e}]  UB={float(UB(k)):+.6e}  Q={float(Qp(1,k)):.5f}  (1/2)P1={float(P(1,k)/2):.6f}")
first = next(k for k in range(2, 400) if UB(k) < 0)
print("\nsmallest k with UB(k)<0:", first, " UB(91)=", UB(91), " UB(92)=", UB(92))
# k <= 91: lower bound G_main - E_bound > 0 ?
lows = [(k, G_main(k) - E_bound(k)) for k in range(2, 92)]
print("min over 2<=k<=91 of lower bound on G(p_k#):", min(lows, key=lambda t: t[1]))
print("all positive:", all(v > 0 for k, v in lows))
# monotonicity sanity of UB
print("UB strictly decreasing k=2..300:", all(UB(k+1) < UB(k) for k in range(2, 300)))
# heuristic constant
C = 1/(120*zeta(4)) - 1/(252*zeta(6)) + 1/(240*zeta(8)) + c4
print("limit constant 1/(120 z4) - 1/(252 z6) + 1/(240 z8) + c4 =", C)
print("6/pi^2 =", 6/pi**2, " flip Q solving (0.5)(6/pi^2/Q)(1-Q/6)+C=0 -> Q =", 1/(mpf(1)/6 - C*pi**2/3))
# non-squarefree: n = 2 * p_k#  (q = 2): G = sum_j c_j q^{-s_j} P_{s_j}, |E| <= (1/132) q^-10 prod(1+p^-10)
def G2_main(k, q=2):
    return P(1,k)/(2*q) - P(2,k)/(12*q**2) + P(4,k)/(120*q**4) - P(6,k)/(252*q**6) + P(8,k)/(240*q**8)
def E2_bound(k, q=2): return Qp(10,k)/(132*q**10)
# scan coarse then fine
import bisect
ks = list(range(1000, 12000, 250))
vals = [(k, G2_main(k)) for k in ks]
for k, v in vals:
    if v < 0:
        print("q=2 coarse: first negative near k=", k, "p_k=", primes[k-1], "Q=", Qp(1,k)); break
