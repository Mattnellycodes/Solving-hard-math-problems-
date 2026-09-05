#!/usr/bin/env python3
"""
Part (1): exact computation of OEIS A067857 for n <= 411 directly from the definition

        sum_{k | n} a(k) / k!  =  H_n  (harmonic number),

i.e. a(n) = n! * ( H_n - sum_{k|n, k<n} a(k)/k! ), computed recursively with exact
fractions, and independently via Moebius inversion a(n) = n! * F(n),
F(n) = sum_{d|n} mu(n/d) H_d.  Both must agree, a(n) must be an integer, and the
conjecture (Robert Israel, 2015)

        a(n) < 0   <=>   omega(n) is odd and omega(n) >= 3

is checked for every n <= 411 (omega = number of DISTINCT prime factors).

Everything here is exact integer/rational arithmetic (Python fractions), no floating point.
"""
import json, os, sys
from fractions import Fraction
from math import factorial

N_MAX = 411
HERE = os.path.dirname(os.path.abspath(__file__))


def smallest_prime_factor_table(n):
    spf = list(range(n + 1))
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


SPF = smallest_prime_factor_table(N_MAX)


def factor(n):
    """Return dict prime -> exponent (exact trial division via the spf table)."""
    f = {}
    while n > 1:
        p = SPF[n]
        while n % p == 0:
            n //= p
            f[p] = f.get(p, 0) + 1
    return f


def omega(n):
    return len(factor(n))


def big_omega(n):
    return sum(factor(n).values())


def mobius(n):
    f = factor(n)
    if any(e > 1 for e in f.values()):
        return 0
    return -1 if len(f) % 2 else 1


def divisors(n):
    return [d for d in range(1, n + 1) if n % d == 0]


# exact harmonic numbers H_0..H_N
H = [Fraction(0)]
for i in range(1, N_MAX + 1):
    H.append(H[-1] + Fraction(1, i))

# (i) recursively from the defining relation sum_{k|n} a(k)/k! = H_n
a_def = {}
for n in range(1, N_MAX + 1):
    s = Fraction(0)
    for k in divisors(n):
        if k < n:
            s += Fraction(a_def[k], factorial(k))
    val = factorial(n) * (H[n] - s)
    assert val.denominator == 1, ("a(%d) is not an integer" % n, val)
    a_def[n] = int(val)

# (ii) via Moebius inversion a(n) = n! * sum_{d|n} mu(n/d) H_d
a_mob = {}
F_exact = {}
for n in range(1, N_MAX + 1):
    F = sum((mobius(n // d) * H[d] for d in divisors(n)), Fraction(0))
    F_exact[n] = F
    val = factorial(n) * F
    assert val.denominator == 1
    a_mob[n] = int(val)

assert a_def == a_mob, "recursive definition and Moebius inversion disagree"

# (iii) re-verify the defining relation for every n with the computed integers
for n in range(1, N_MAX + 1):
    assert sum((Fraction(a_def[k], factorial(k)) for k in divisors(n)), Fraction(0)) == H[n]

# anchor values (OEIS A067857 / formal-conjectures test theorems)
anchors = {1: 1, 2: 1, 3: 5, 4: 14, 5: 154, 30: -22690644647302814715858124800000}
for n, v in anchors.items():
    assert a_def[n] == v, (n, a_def[n], v)

# (iv) the conjecture for n <= 411 (omega = distinct prime factors), and the
#      mis-formalised Omega (with multiplicity) variant for comparison
violations_omega = []
violations_Omega = []
negatives = []
for n in range(1, N_MAX + 1):
    neg = a_def[n] < 0
    w = omega(n)
    W = big_omega(n)
    if neg:
        negatives.append(n)
    if neg != (w % 2 == 1 and w >= 3):
        violations_omega.append(n)
    if neg != (W % 2 == 1 and W >= 3):
        violations_Omega.append(n)

print("A067857 exact, n <= %d" % N_MAX)
print("a(1..12) =", [a_def[n] for n in range(1, 13)])
print("a(30) =", a_def[30])
print("first negative index:", negatives[0], "; all n<=411 with a(n)<0:", negatives)
print("omega parity check: violations of 'a(n)<0 <=> omega(n) odd and >=3' for n<=411:", violations_omega)
print("(for comparison) violations of the Omega-with-multiplicity variant for n<=411:", violations_Omega[:20],
      "... (%d in total, first one n=%s)" % (len(violations_Omega), violations_Omega[0] if violations_Omega else None))
print("check: every n<=411 with omega(n) odd >=3 has a(n)<0:",
      all(a_def[n] < 0 for n in range(1, N_MAX + 1) if omega(n) % 2 == 1 and omega(n) >= 3))
print("check: every n<=411 with omega(n) even or omega(n)=1 has a(n)>0:",
      all(a_def[n] > 0 for n in range(1, N_MAX + 1) if not (omega(n) % 2 == 1 and omega(n) >= 3)))
print("a(n) == 0 for some n<=411?", any(a_def[n] == 0 for n in range(1, N_MAX + 1)))

# sign of F(p_k#) for the small primorials that are <= 411: 2, 6, 30, 210
for n in (2, 6, 30, 210):
    print("n=%d omega=%d F(n)=%s  a(n)=%d" % (n, omega(n), F_exact[n], a_def[n]))

out = {
    "sequence": "A067857",
    "definition": "sum_{k|n} a(k)/k! = H_n ; equivalently a(n) = n! * sum_{d|n} mu(n/d) H_d",
    "n_max": N_MAX,
    "arithmetic": "exact rational (python fractions); both definitions agree and satisfy the defining relation",
    "a_first_40": {str(n): a_def[n] for n in range(1, 41)},
    "a_30": a_def[30],
    "negative_indices_le_411": negatives,
    "conjecture_omega_violations_le_411": violations_omega,
    "conjecture_Omega_with_multiplicity_violations_le_411": violations_Omega,
    "F_exact_primorials": {str(n): [F_exact[n].numerator, F_exact[n].denominator] for n in (2, 6, 30, 210)},
}
with open(os.path.join(HERE, "results_part1_small_n.json"), "w") as f:
    json.dump(out, f, indent=1)
print("wrote results_part1_small_n.json")
