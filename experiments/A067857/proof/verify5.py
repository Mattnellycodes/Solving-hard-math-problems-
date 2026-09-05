# incremental (O(k)) scan: G(2*p_k#) lower bound positive for all 2<=k<=9230; and G(p_k#) enclosure summary; digits.
from sympy import primerange
from mpmath import mp, mpf, euler, pi, harmonic, log, bernoulli
import math
mp.dps = 30
primes = list(primerange(2, 100000))
R42 = harmonic(2) - log(2) - euler - (mpf(1)/4 - sum(bernoulli(2*j)/(2*j*mpf(2)**(2*j)) for j in range(1,5)))
Ecap2 = mpf(1)/(132*1024)*mpf('1e-3')
P = {s: mpf(1) for s in (1,2,4,6,8)}
minval = None; first_neg = None
for k in range(1, 9400):
    p = primes[k-1]
    for s in P: P[s] *= (1 - mpf(p)**(-s))
    if k < 2: continue
    g = P[1]/4 - P[2]/48 + P[4]/1920 - P[6]/16128 + P[8]/61440 + R42
    if k <= 9230:
        lb = g - Ecap2
        if minval is None or lb < minval[1]: minval = (k, lb)
    elif first_neg is None and g + Ecap2 < 0: first_neg = (k, g + Ecap2)
print("n_k = 4*3*5*...*p_k: min over 2<=k<=9230 of certified lower bound on G(n_k):", minval, " -> all positive:", minval[1] > 0)
print("first k with certified G(n_k)<0:", first_neg, " p_k =", primes[first_neg[0]-1])
print("log10(2 * p_9231#) =", sum(math.log10(q) for q in primes[:9231]) + math.log10(2))
print("log10(479#) =", sum(math.log10(q) for q in primes[:92]), " log10(487#)=", sum(math.log10(q) for q in primes[:93]))
# squarefree primorials: incremental enclosure summary for k=2..100 (m=4 formula, per-n error bound)
P = {s: mpf(1) for s in (1,2,4,6,8,10)}; Q10 = mpf(1)
c4 = mpf(2897)/5040 - euler
rows = []
for k in range(1, 101):
    p = primes[k-1]
    for s in P: P[s] *= (1 - mpf(p)**(-s))
    Q10 *= (1 + mpf(p)**(-10))
    if k < 2: continue
    g = P[1]/2 - P[2]/12 + P[4]/120 - P[6]/252 + P[8]/240 + c4; e = (Q10 - 1)/132
    rows.append((k, g - e, g + e))
print("G(p_k#) enclosures: min lower bound for k<=91:", min((r for r in rows if r[0] <= 91), key=lambda r: r[1]))
print("k=92:", rows[90], " k=93:", rows[91])
print("all k in 2..91 certified positive:", all(r[1] > 0 for r in rows if r[0] <= 91), "; all k in 92..100 certified negative:", all(r[2] < 0 for r in rows if r[0] >= 92))
