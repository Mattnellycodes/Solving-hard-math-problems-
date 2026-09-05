#!/usr/bin/env python3
"""Independent exact recomputation of A067857 for n <= 411 (auditor's code, sympy number theory + Fractions),
compared with (i) the 20 OEIS data terms (from github.com/oeis/oeisdata seq/A067/A067857.seq, fetched today),
(ii) the certificate's results_part1_small_n.json, (iii) the Lean definition a(n) = n! sum_{d|n} mu(n/d) H_d."""
import json, sys
from fractions import Fraction
from math import factorial
from sympy import divisors, mobius, primefactors, factorint, harmonic, Rational
sys.set_int_max_str_digits(0)
NMAX = 411
OEIS = [1,1,5,14,154,84,8028,25584,361296,528480,80627040,33471360,13575738240,13835646720,263577888000,
        13869128448000,867718162483200,316745643110400,309920046408806400,207862451693568000]
H = [Fraction(0)]
for i in range(1, NMAX + 1):
    H.append(H[-1] + Fraction(1, i))
# sympy's harmonic as an independent source for H_d
for d in (1, 2, 30, 100, 411):
    hs = harmonic(d); assert Fraction(int(hs.p), int(hs.q)) == H[d]
a = {}
for n in range(1, NMAX + 1):
    F = sum((Fraction(int(mobius(n // d))) * H[d] for d in divisors(n)), Fraction(0))
    v = factorial(n) * F
    assert v.denominator == 1
    a[n] = int(v)
# defining relation sum_{k|n} a(k)/k! = H_n
for n in range(1, NMAX + 1):
    assert sum((Fraction(a[k], factorial(k)) for k in divisors(n)), Fraction(0)) == H[n]
assert [a[n] for n in range(1, 21)] == OEIS, "mismatch with OEIS data terms"
assert a[30] == -22690644647302814715858124800000
J = json.load(open("/home/user/Solving-hard-math-problems-/experiments/A067857/certificate/results_part1_small_n.json"))
assert all(a[int(k)] == v for k, v in J["a_first_40"].items())
neg = [n for n in range(1, NMAX + 1) if a[n] < 0]
assert neg == J["negative_indices_le_411"], "negative index list differs from certificate"
viol = [n for n in range(1, NMAX + 1) if (a[n] < 0) != (len(primefactors(n)) % 2 == 1 and len(primefactors(n)) >= 3)]
violOmega = [n for n in range(1, NMAX + 1) if (a[n] < 0) != (sum(factorint(n).values()) % 2 == 1 and sum(factorint(n).values()) >= 3)]
print("a(1..20) == OEIS data terms: True;  a(30) == OEIS value: True")
print("certificate a_first_40 and negative-index list (%d entries) reproduced exactly" % len(neg))
print("violations of 'a(n)<0 <=> omega(n) odd >= 3' for n <= 411:", viol)
print("violations of the Omega (with multiplicity) variant, first 10:", violOmega[:10], "count", len(violOmega))
print("any a(n) == 0 for n <= 411:", any(a[n] == 0 for n in a))
print("F(30) = a(30)/30! =", Fraction(a[30], factorial(30)), " F(210) sign:", 1 if a[210] > 0 else -1, " omega(210) = 4")
json.dump({"n_max": NMAX, "negative_indices": neg, "violations_omega": viol, "violations_Omega": violOmega,
           "a_first_20_match_oeis": True, "matches_certificate_json": True}, open("results_small_n.json", "w"), indent=1)
