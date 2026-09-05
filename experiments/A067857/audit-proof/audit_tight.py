#!/usr/bin/env python3
# Independent tight rigorous enclosure of G(p_k#) (k = 2..12 and 88..96) by a different route than PROOF.md's coarse bound:
# closed form (m=4) + exact per-divisor corrections eps(d)-A_4(d) for all divisors d<=D (interval H_d, log d, gamma via MPFR
# directed rounding) + first-omitted-term tail bound for divisors d>D. Compared with exact rationals for k<=7.
import time
from ivmpfr import *
from sympy import primerange
t0 = time.time()
D = 60000
primes = list(primerange(2, 1300))
Hlo = [mpfr(0)]; Hhi = [mpfr(0)]
for d in range(1, D+1):
    Hlo.append(dn(lambda: Hlo[-1] + mpfr(1)/d)); Hhi.append(up(lambda: Hhi[-1] + mpfr(1)/d))
def eps(d): return Iv(Hlo[d], Hhi[d]) - LOG(d) - GAMMA
def A4(d): return Fraction(1, 2*d) - Fraction(1, 12*d**2) + Fraction(1, 120*d**4) - Fraction(1, 252*d**6) + Fraction(1, 240*d**8)
corr_cache = {}
def corr_term(d):   # eps(d) - A4(d) as interval (cached)
    if d not in corr_cache: corr_cache[d] = eps(d) - Iv.of(A4(d))
    return corr_cache[d]
print("H_D width:", float(up(lambda: Hhi[D]-Hlo[D])), " [%.1fs]" % (time.time()-t0))
def G_enclosure(k):
    ps = primes[:k]
    P = {s: Iv.of(1) for s in (1,2,4,6,8)}; Q10 = Iv.of(1)
    for p in ps:
        for s in P: P[s] = P[s] * (1 - Iv.of(Fraction(1, p**s)))
        Q10 = Q10 * (1 + Iv.of(Fraction(1, p**10)))
    closed = P[1]/2 - P[2]/12 + P[4]/120 - P[6]/252 + P[8]/240
    # divisors d<=D of p_k#: DFS
    divs = [(1, 0)]
    for p in ps:
        divs += [(d*p, c+1) for d, c in divs if d*p <= D]
    corr = Iv.of(0); S10 = Fraction(0)
    for d, c in divs:
        t = corr_term(d); corr = corr + (t if c % 2 == 0 else -t); S10 += Fraction(1, d**10)
    tail = ((Q10 - Iv.of(S10)) * Iv.of(Fraction(1, 132))).hi
    if tail < 0: tail = mpfr(0)
    tot = closed + corr
    return Iv(dn(lambda: tot.lo - tail), up(lambda: tot.hi + tail)), len(divs), tail
# exact rationals for k<=7
from sympy import divisors, mobius
def Fexact(n):
    Hn = [Fraction(0)]
    for i in range(1, n+1): Hn.append(Hn[-1] + Fraction(1, i))
    return sum(int(mobius(n//d)) * Hn[d] for d in divisors(n))
for k in range(2, 8):
    n = 1
    for p in primes[:k]: n *= p
    G, nd, tail = G_enclosure(k)
    if n <= 30030:
        Ge = (-1)**k * Fexact(n); inside = (G.lo <= mpfr(Ge.numerator)/mpfr(Ge.denominator) <= G.hi)
        print(f"k={k} n={n}: G in {G} width={float(G.width()):.1e}  exact G={float(Ge):+.15e} inside: {inside}")
    else:
        print(f"k={k} n={n}: G in {G} width={float(G.width()):.1e}")
for k in list(range(88, 97)) + [120, 200]:
    G, nd, tail = G_enclosure(k)
    print(f"k={k} p_k={primes[k-1]}: G(p_k#) in {G}  width={float(G.width()):.1e}  #divisors<=D: {nd}  tail={float(tail):.1e}  sign={'+' if G.lo>0 else ('-' if G.hi<0 else '?')}  [%.1fs]" % (time.time()-t0))
