"""Main agent's own sanity check of the A067857 claim (2026-09-05).
a(n) = n! * F(n), F(n) = sum_{d|n} mu(n/d) H_d.  Conjecture: a(n) < 0 iff omega(n) odd and >= 3.
Part 1: exact check for n <= 300.  Part 2: G(k) = (-1)^k F(p_k#) via the closed form with mpmath at 60 digits
(NOT a certificate, just a sanity check of the sign flip near k = 92), and a brute-force comparison for k <= 12.
"""
from fractions import Fraction
from sympy import factorint, primerange, divisors, mobius, primorial
import mpmath as mp, itertools, math
# Part 1
H = [Fraction(0)]
for i in range(1, 301): H.append(H[-1] + Fraction(1, i))
bad = []
for n in range(1, 301):
    F = sum(mobius(n // d) * H[d] for d in divisors(n))
    om = len(factorint(n))
    neg = F < 0
    pred = (om % 2 == 1 and om >= 3)
    if neg != pred: bad.append((n, om, float(F)))
print("n<=300: conjecture violations:", bad)
n30 = 30; F30 = sum(mobius(30 // d) * H[d] for d in divisors(30)); print("a(30) =", math.factorial(30) * F30)
# Part 2
mp.mp.dps = 60
def eps(d):  # H_d - ln d - gamma exactly via digamma: H_d = psi(d+1) + gamma
    return mp.psi(0, d + 1) - mp.log(d)
primes = list(primerange(2, 700))
def G_closed(k, m=6, D=20000):
    ps = primes[:k]
    # eps(d) = 1/(2d) - sum_{j=1..m} B_{2j}/(2j d^{2j}) + R_m(d)
    # sum_{d|n} mu(n/d) d^{-s} = (-1)^k prod (1 - p^{-s})
    coef = {1: mp.mpf(1) / 2}
    for j in range(1, m + 1):
        coef[2 * j] = -mp.bernoulli(2 * j) / (2 * j)
    total = mp.mpf(0)
    for s, c in coef.items():
        prod = mp.mpf(1)
        for p in ps: prod *= (1 - mp.mpf(p) ** (-s))
        total += c * prod
    # exact correction for divisors d <= D: add mu-weighted (eps(d) - A_m(d)), sign (-1)^k folded: G = (-1)^k F
    # F = sum_d mu(n/d) eps(d); with n squarefree, mu(n/d) = (-1)^{k - omega(d)}; G = (-1)^k F = sum_d (-1)^{omega(d)} eps(d)
    # so G = sum_d (-1)^{omega(d)} [A_m(d) + R_m(d)] = closed-form part + sum_d (-1)^{omega(d)} R_m(d)
    def A(d):
        return sum(c * mp.mpf(d) ** (-s) for s, c in coef.items())
    corr = mp.mpf(0); count = 0
    # enumerate squarefree divisors d <= D of p_k#
    def rec(i, d, om):
        nonlocal corr, count
        corr += (-1) ** om * (eps(d) - A(d)); count += 1
        for j in range(i, k):
            if d * ps[j] > D: break
            rec(j + 1, d * ps[j], om + 1)
    rec(0, 1, 0)
    tail = abs(mp.bernoulli(2 * m + 2)) / (2 * m + 2) * mp.mpf(D) ** (-(2 * m + 2)) * (2 ** k)  # crude bound on ignored remainders
    return total + corr, tail, count
def G_brute(k):
    ps = primes[:k]; tot = mp.mpf(0)
    for r in range(k + 1):
        for sub in itertools.combinations(ps, r):
            d = 1
            for p in sub: d *= p
            tot += (-1) ** r * eps(d)
    return tot
for k in (6, 9, 12):
    gc, tail, cnt = G_closed(k, m=4, D=10**6)
    print(f"k={k}: closed={mp.nstr(gc, 25)} brute={mp.nstr(G_brute(k), 25)} diff={mp.nstr(gc - G_brute(k), 5)}")
for k in (88, 90, 91, 92, 93, 94, 95):
    gc, tail, cnt = G_closed(k, m=6, D=20000)
    print(f"k={k}: G ≈ {mp.nstr(gc, 12)}  (crude tail bound {mp.nstr(tail, 3)}; {cnt} divisors corrected)  => a(p_k#) sign = {'+' if (gc > 0) == (k % 2 == 0) else '-'}  conjecture predicts {'-' if k % 2 == 1 else '+'}")
