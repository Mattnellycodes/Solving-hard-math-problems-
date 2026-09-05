# Independent verification, part 1: exact-rational F(n) vs. the closed form (m=4), small squarefree n.
from fractions import Fraction
from math import factorial
from sympy import primerange, divisors, mobius, primefactors, bernoulli as sbern
from mpmath import mp, mpf, euler, log, harmonic, fprod
mp.dps = 40

def F_exact(n):
    H = [Fraction(0)]
    for i in range(1, n+1): H.append(H[-1] + Fraction(1, i))
    return sum(int(mobius(n//d)) * H[d] for d in divisors(n))

# OEIS cross-check: a(30) and small-n sign pattern
print("a(30) =", factorial(30)*F_exact(30), " (OEIS: -22690644647302814715858124800000)")
viol = [n for n in range(1, 200) if (F_exact(n) < 0) != (len(primefactors(n)) % 2 == 1 and len(primefactors(n)) >= 3)]
print("n<200 violating Israel's conjecture:", viol)

# closed form for squarefree n with prime set ps, m=4:
# G(n) = (-1)^omega F(n) = 1/2 P1 - 1/12 P2 + 1/120 P4 - 1/252 P6 + 1/240 P8 + c4 + E4,  |E4| <= (1/132)(prod(1+p^-10)-1)
c4 = Fraction(2897, 5040)  # 1/2+1/12-1/120+1/252-1/240
assert c4 == Fraction(1,2)+Fraction(1,12)-Fraction(1,120)+Fraction(1,252)-Fraction(1,240)
def P(ps, s): return fprod([1 - mpf(p)**(-s) for p in ps])
def G_main(ps):
    return P(ps,1)/2 - P(ps,2)/12 + P(ps,4)/120 - P(ps,6)/252 + P(ps,8)/240 + (mpf(c4.numerator)/c4.denominator - euler)
def E_bound(ps): return (fprod([1 + mpf(p)**(-10) for p in ps]) - 1)/132
import itertools
primes = list(primerange(2, 60))
worst = 0
for k in range(2, 6):
    for ps in itertools.combinations(primes[:9], k):
        n = 1
        for p in ps: n *= p
        G_ex = (-1)**k * F_exact(n)
        G_ex = mpf(G_ex.numerator)/G_ex.denominator
        diff = abs(G_ex - G_main(ps)); eb = E_bound(ps)
        worst = max(worst, diff/eb)
        assert diff <= eb, (ps, diff, eb)
print("all squarefree n with 2..5 primes from first 9 primes: |G_exact - G_main| <= E_bound; max ratio |diff|/E_bound =", worst)
# primorials k=2..8 exact
for k in range(2, 9):
    ps = primes[:k]; n = 1
    for p in ps: n *= p
    G_ex = (-1)**k * F_exact(n); G_ex = mpf(G_ex.numerator)/G_ex.denominator
    print(f"k={k} n={n}: G_exact={G_ex}  G_main={G_main(ps)}  diff={G_ex-G_main(ps)}  E_bound={E_bound(ps)}")
