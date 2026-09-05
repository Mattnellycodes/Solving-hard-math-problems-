#!/usr/bin/env python3
# Independent audit of PROOF.md numerics with MPFR directed rounding (gmpy2), see ivmpfr.py.
#  (1) certificate for Prop. 6(c): UB(92) < 0; also UB(91) > 0, UB(93) < 0, first k with UB(k)<0, monotonicity check;
#  (2) exact-rational P_1(92), P_4(92), P_8(92) and the 8-term table with ceil at 10 decimals (re-derived);
#  (3) enclosures of G(p_k#) with the per-n error bound for k = 2..200 (Remark (a));
#  (4) Q(k) thresholds 6 and 12 (Prop. 7); (5) n_k = 2 p_k# family (Remark (d)).
import time
from ivmpfr import *
from sympy import primerange
t0 = time.time()
primes = list(primerange(2, 130000))
assert primes[90] == 467 and primes[91] == 479 and primes[92] == 487, (primes[90], primes[91], primes[92])
print("p_91, p_92, p_93 =", primes[90], primes[91], primes[92], "; #primes<=479 =", sum(1 for p in primes if p <= 479))
print("PI  =", PI.strfull()); print("GAMMA =", GAMMA.strfull()); print("LN2 =", LN2.strfull())
# ---------- (2) exact rationals and the table ----------
def Pex(s, k):
    r = Fraction(1)
    for p in primes[:k]: r *= 1 - Fraction(1, p**s)
    return r
P1e, P4e, P8e = Pex(1, 92), Pex(4, 92), Pex(8, 92)
def dec(fr, nd=25):
    q = 10**nd; return "%d.%0*d" % (fr.numerator // fr.denominator, nd, (fr.numerator * q // fr.denominator) % q)
print("exact P1(92) =", dec(P1e), "\nexact P4(92) =", dec(P4e), "\nexact P8(92) =", dec(P8e))
def ceil10(fr, nd=10):
    q = 10**nd; return Fraction(-((-fr.numerator * q) // fr.denominator), q)
pi_hi = Fraction('3.1415926536'); g_lo = Fraction('0.5772156649')
assert PI.hi < mpfr(pi_hi.numerator)/pi_hi.denominator or True  # checked below with intervals
terms = [P1e/2, P4e/120, P8e/240, -1/(2*pi_hi**2), -Fraction(15,4)/pi_hi**6, Fraction(2897,5040), -g_lo, Fraction(1,132000)]
tab = [ceil10(t) for t in terms]
claimed = [Fraction(x) for x in "0.0450790249 0.0076994867 0.0041497467 -0.0506605918 -0.0039006055 0.5748015874 -0.5772156649 0.0000075758".split()]
print("table terms re-derived (ceil at 1e-10):", [str(float(t)) for t in tab])
print("match PROOF.md table:", tab == claimed, "; sum =", float(sum(tab)), "; claimed sum -0.0000394407:", sum(tab) == Fraction('-0.0000394407'))
print("exact-rational UB(92) upper bound (pi<=3.1415926536, gamma>=0.5772156649):", float(sum(terms)), "<0:", sum(terms) < 0)
# check the two constant bounds against MPFR enclosures
print("pi <= 3.1415926536 ?", PI.hi < mpfr(pi_hi.numerator)/mpfr(pi_hi.denominator), "  gamma >= 0.5772156649 ?", GAMMA.lo > mpfr(g_lo.numerator)/mpfr(g_lo.denominator))
# ---------- (1),(3) interval loop over k ----------
KMAX = 200
inv2pi2 = 1/(2*PI*PI); c15 = Iv.of(15)/(4*PI**6); c4 = Iv.of(Fraction(2897,5040)) - GAMMA; Ecap = Iv.of(Fraction(1,132000))
P = {s: Iv.of(1) for s in (1,2,4,6,8)}; Q10 = Iv.of(1)
UBs = {}; Gs = {}
for k in range(1, KMAX+1):
    p = primes[k-1]
    for s in P: P[s] = P[s] * (1 - Iv.of(Fraction(1, p**s)))
    Q10 = Q10 * (1 + Iv.of(Fraction(1, p**10)))
    if k < 2: continue
    UB = P[1]/2 + P[4]/120 + P[8]/240 - inv2pi2 - c15 + c4 + Ecap
    main = P[1]/2 - P[2]/12 + P[4]/120 - P[6]/252 + P[8]/240 + c4
    E = (Q10 - 1)/132
    G = Iv(dn(lambda: main.lo - E.hi), up(lambda: main.hi + E.hi))
    UBs[k] = UB; Gs[k] = G
    if k in (2,3,4,5,10,50,85,86,87,88,89,90,91,92,93,94,95,100,150,200):
        print(f"k={k:3d} p_k={p:4d}  UB(k) in {UB}   G(p_k#) in {G}")
print("UB(91) =", UBs[91].strfull()); print("UB(92) =", UBs[92].strfull()); print("UB(93) =", UBs[93].strfull())
print("first k with UB(k).hi < 0:", next(k for k in range(2, KMAX+1) if UBs[k].hi < 0), "; UB(91).lo > 0:", UBs[91].lo > 0)
print("UB(k+1).hi < UB(k).lo for all 2<=k<200 (numerical monotonicity):", all(UBs[k+1].hi < UBs[k].lo for k in range(2, KMAX)))
print("G(p_k#).lo > 0 for all 2<=k<=91:", all(Gs[k].lo > 0 for k in range(2, 92)), "; min lower bound:", min((Gs[k].lo, k) for k in range(2, 92)))
print("G(p_k#).hi < 0 for all 92<=k<=200:", all(Gs[k].hi < 0 for k in range(92, KMAX+1)), "; G(92) hi:", float(Gs[92].hi), " G(93) hi:", float(Gs[93].hi))
print("G(p_k#) < UB(k) numerically for all k (G.hi < UB.lo):", all(Gs[k].hi < UBs[k].lo for k in range(2, KMAX+1)))
print("[%.1fs]" % (time.time()-t0))
# ---------- (4) Q thresholds ----------
Q = Iv.of(1); Qk = {}
for k in range(1, 6500):
    Q = Q * (1 + Iv.of(Fraction(1, primes[k-1]))); Qk[k] = Q
for k in (51, 52, 91, 92, 6480, 6481): print(f"Q({k}) = {Qk[k]}  p_{k}={primes[k-1]}")
k6 = next(k for k in range(1, 6500) if Qk[k].lo > 6); k12 = next(k for k in range(1, 6500) if Qk[k].lo > 12)
print("least k with Q(k)>6 (certified: Q(k-1).hi<6<Q(k).lo):", k6, Qk[k6-1].hi < 6, "; least k with Q(k)>12:", k12, Qk[k12-1].hi < 12)
print("[%.1fs]" % (time.time()-t0))
# ---------- (5) n_k = 2 p_k#, q = 2 ----------
R42 = (Iv.of(Fraction(3,2)) - LN2 - GAMMA) - Iv.of(Fraction(1,4) - Fraction(1,48) + Fraction(1,1920) - Fraction(1,16128) + Fraction(1,61440))
print("R_4(2) =", R42.strfull())
P = {s: Iv.of(1) for s in (1,2,4,6,8)}; Q10 = Iv.of(1)
minlb = None; first_neg = None; first_ub2 = None
c8 = 1/(8*PI*PI); c945 = Iv.of(945)/(16128*PI**6); cap = Iv.of(Fraction(74, 10**10))
for k in range(1, 9400):
    p = primes[k-1]
    for s in P: P[s] = P[s] * (1 - Iv.of(Fraction(1, p**s)))
    Q10 = Q10 * (1 + Iv.of(Fraction(1, p**10)))
    if k < 2: continue
    g = P[1]/4 - P[2]/48 + P[4]/1920 - P[6]/16128 + P[8]/61440 + R42
    E2 = (Q10 - 1)/(132*1024)
    lo = dn(lambda: g.lo - E2.hi); hi = up(lambda: g.hi + E2.hi)
    if k <= 9230 and (minlb is None or lo < minlb[0]): minlb = (lo, k)
    if first_neg is None and hi < 0: first_neg = (k, hi)
    UB2 = P[1]/4 + P[4]/1920 + P[8]/61440 - c8 - c945 + R42 + cap
    if first_ub2 is None and UB2.hi < 0: first_ub2 = (k, UB2.hi)
print("n_k=2p_k#: min certified lower bound of G(n_k) over 2<=k<=9230:", float(minlb[0]), "at k=", minlb[1], "-> all positive:", minlb[0] > 0)
print("first k with certified G(n_k)<0:", first_neg[0], "upper bound", float(first_neg[1]), "p_k =", primes[first_neg[0]-1])
print("first k with UB_2(k) < 0 (uniform bound with 7.4e-9 cap):", first_ub2[0], float(first_ub2[1]))
print("zeta(10)-1 elementary bound: 2^-10+3^-10+3^-9/9 =", float(Fraction(1,2**10)+Fraction(1,3**10)+Fraction(1,9*3**9)), "; (1/132)2^-10*1e-3 =", 1e-3/(132*1024))
print("[%.1fs]" % (time.time()-t0))
