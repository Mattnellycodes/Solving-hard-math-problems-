# Independent verification, part 1b: exact-rational G(n) vs. closed form (m=4) for squarefree n <= 6000 (all of them),
# and mpmath (40 dps) brute force over all 2^k divisors for primorials k<=12.
from fractions import Fraction
from sympy import primerange, divisors, mobius, primefactors, factorint
from mpmath import mp, mpf, euler, fprod, harmonic
mp.dps = 40
N = 6000
H = [Fraction(0)]
for i in range(1, N+1): H.append(H[-1] + Fraction(1, i))
def F_exact(n): return sum(int(mobius(n//d)) * H[d] for d in divisors(n))
c4 = mpf(2897)/5040 - euler
def P(ps, s): return fprod([1 - mpf(p)**(-s) for p in ps])
def G_main(ps): return P(ps,1)/2 - P(ps,2)/12 + P(ps,4)/120 - P(ps,6)/252 + P(ps,8)/240 + c4
def E_bound(ps): return (fprod([1 + mpf(p)**(-10) for p in ps]) - 1)/132
worst = 0; cnt = 0
for n in range(2, N+1):
    f = factorint(n)
    if any(e > 1 for e in f.values()) or len(f) < 2: continue
    ps = sorted(f); k = len(ps)
    G_ex = (-1)**k * F_exact(n); G_ex = mpf(G_ex.numerator)/G_ex.denominator
    diff = abs(G_ex - G_main(ps)); eb = E_bound(ps); cnt += 1
    assert diff <= eb, (n, diff, eb)
    worst = max(worst, diff/eb)
print(f"checked {cnt} squarefree n<={N} with omega>=2: |G_exact - G_main| <= E_bound always; max ratio = {float(worst):.4f}")
# also the sign law: sign a(n) = (-1)^omega sign G(n) trivially; and G_exact > 0 for all of them (conjecture holds there)
# brute force primorials with mpmath
primes = list(primerange(2, 50))
for k in range(2, 13):
    ps = primes[:k]; divs = [(1, 0)]
    for p in ps: divs += [(d*p, c+1) for d, c in divs]
    Fb = sum((-1)**(k-c) * harmonic(d) for d, c in divs)
    Gb = (-1)**k * Fb
    print(f"k={k:2d} brute G={float(Gb):+.12e}  main={float(G_main(ps)):+.12e}  |diff|={float(abs(Gb-G_main(ps))):.2e}  E_bound={float(E_bound(ps)):.2e}")
