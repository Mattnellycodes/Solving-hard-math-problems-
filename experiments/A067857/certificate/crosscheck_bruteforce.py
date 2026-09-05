#!/usr/bin/env python3
"""
Cross-check of the closed-form certificates against brute force over ALL 2^k divisors of p_k#.

  (i)  k <= 7  : F(p_k#) = sum_{d|n} mu(n/d) H_d computed EXACTLY (rational), via an lcm
                 accumulator for H_d, d <= n = p_k#  (n = 510510 for k = 7).  The exact value
                 must lie INSIDE both certificate enclosures (rigorous containment test).
  (ii) k <= 14 : brute force with mpmath at 60 digits (mp.harmonic), 2^14 = 16384 divisors.
                 Reported: distance of the brute-force value from each enclosure (0 = inside);
                 this is a floating-point consistency check, not a rigorous one.
"""
import json, os, sys, time
from fractions import Fraction
from mpmath import mp, mpf, harmonic, fsum

sys.set_int_max_str_digits(0)
HERE = os.path.dirname(os.path.abspath(__file__))
mp.dps = 60

def first_primes(k):
    ps, q = [], 2
    while len(ps) < k:
        if all(q % p for p in ps):
            ps.append(q)
        q += 1
    return ps

def all_divisors(ps):
    divs = [(1, 0)]
    for p in ps:
        divs += [(d * p, w + 1) for d, w in divs]
    return divs

def exact_F_primorial(k):
    """Exact F(n), n = p_k#, via H_d = N_d / L_d, L_d = lcm(1..d), accumulated without gcds."""
    ps = first_primes(k)
    n = 1
    for p in ps:
        n *= p
    signs = {d: (-1) ** (k - w) for d, w in all_divisors(ps)}        # mu(n/d) = (-1)^{omega(n/d)} = (-1)^{k - omega(d)}
    spf = list(range(n + 1))
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    N, L, S = 0, 1, 0                                                  # H_d = N/L ; S/L = running sum over divisors
    for d in range(1, n + 1):
        if d > 1:
            q, x = spf[d], d
            while x % q == 0:
                x //= q
            if x == 1:                                                 # d is a prime power -> lcm grows by q
                L *= q; N *= q; S *= q
        N += L // d
        sg = signs.get(d)
        if sg is not None:
            S += sg * N
    return Fraction(S, L)

certs = {}
for route in ("a", "b"):
    J = json.load(open(os.path.join(HERE, "cert_route_%s.json" % route)))
    certs[route] = {c["k"]: (Fraction(c["G_enclosure"][0]), Fraction(c["G_enclosure"][1])) for c in J["certificates"]}

results = {"exact": {}, "mpmath60": {}}
print("(i) exact rational brute force, k <= 7 (containment in the certified enclosures is required)")
for k in range(2, 8):
    t = time.time()
    Fex = exact_F_primorial(k)
    Gex = (-1) ** k * Fex
    line = "k=%d n=%-7d G exact = %+.18e" % (k, eval("*".join(map(str, first_primes(k)))), float(Gex))
    ok_all = True
    for route in ("a", "b"):
        lo, hi = certs[route][k]
        inside = lo <= Gex <= hi
        ok_all &= inside
        line += "  route %s: %s (margin %.1e)" % (route, "inside" if inside else "OUTSIDE", float(min(Gex - lo, hi - Gex)))
    print(line + "   [%.1fs]" % (time.time() - t))
    results["exact"][k] = {"G_exact": [Gex.numerator, Gex.denominator] if k <= 4 else "%.40e" % float(Gex),
                           "inside_route_a": certs["a"][k][0] <= Gex <= certs["a"][k][1],
                           "inside_route_b": certs["b"][k][0] <= Gex <= certs["b"][k][1]}
    assert ok_all, "exact brute force lies outside a certified enclosure at k=%d" % k
    if k == 2:
        assert Fex == Fraction(7, 60)              # matches part (1): F(6) = 7/60
    if k == 3:
        assert Fex == Fraction(-199238355653, 2329089562800)   # matches part (1): F(30)

print("\n(ii) 60-digit mpmath brute force over all 2^k divisors, k <= 14 (distance 0 = inside enclosure)")
for k in range(2, 15):
    t = time.time()
    ps = first_primes(k)
    G = fsum([(-1) ** w * harmonic(d) for d, w in all_divisors(ps)])
    Gf = Fraction(mp.nstr(G, 58, strip_zeros=False))   # exact decimal parse of the 58-digit value
    line = "k=%2d  G (brute, 60 dps) = %s" % (k, mp.nstr(G, 25))
    dists = {}
    for route in ("a", "b"):
        lo, hi = certs[route][k]
        dist = max(Fraction(0), lo - Gf, Gf - hi)
        dists[route] = float(dist)
        line += "   dist to route %s enclosure: %.1e" % (route, float(dist))
    print(line + "   [%.1fs]" % (time.time() - t))
    results["mpmath60"][k] = {"G_brute": mp.nstr(G, 45), "dist_route_a": dists["a"], "dist_route_b": dists["b"]}
    assert max(dists.values()) < 1e-45, "brute force disagrees with certificate at k=%d" % k

with open(os.path.join(HERE, "results_crosscheck_bruteforce.json"), "w") as f:
    json.dump(results, f, indent=1)
print("\nall brute-force cross-checks consistent; wrote results_crosscheck_bruteforce.json")
