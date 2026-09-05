#!/usr/bin/env python3
"""Cross-examination of the certificate JSON files themselves (cert_route_a.json, cert_route_b.json) with my own code:
 (1) per-k divisor counts 'num_small_divisors_used' vs an independent count (spf sieve over all integers <= D);
 (2) every stored Euler-product enclosure must contain the exact rational P_s(k);
 (3) their tail bounds T(k) recomputed EXACTLY with their own (m, D): route b stores the exact tail sum, route a an upper bound;
 (4) sign/parity/status bookkeeping for k = 91, 92, 93 and margin |G| / width;
 (5) the JSON G-enclosures must contain my route-C enclosures (already asserted in audit_routeC.py; repeated here from the files)."""
import json, sys
from fractions import Fraction
from math import comb
sys.set_int_max_str_digits(0)
CERT = "/home/user/Solving-hard-math-problems-/experiments/A067857/certificate"
HERE = "/home/user/Solving-hard-math-problems-/experiments/A067857/audit-certificate"
mine = {r["k"]: r for r in json.load(open(HERE + "/results_routeC.json"))["results"]}
def bern(nmax):
    Bs = [Fraction(1)]
    for n in range(1, nmax + 1):
        Bs.append(-sum(comb(n + 1, j) * Bs[j] for j in range(n)) / (n + 1))
    return Bs
for route in ("a", "b"):
    J = json.load(open("%s/cert_route_%s.json" % (CERT, route)))
    m, D, KMAX = J["parameters"]["m"], J["parameters"]["D"], J["parameters"]["KMAX"]
    primes = J["primes"]; pidx = {p: i + 1 for i, p in enumerate(primes)}
    S_TAIL = 2 * m + 2
    Bs = bern(S_TAIL); C_TAIL = abs(Bs[S_TAIL]) / S_TAIL
    assert Fraction(J["tail_constant"]) == C_TAIL
    # (1) independent divisor enumeration: spf sieve
    spf = list(range(D + 1))
    for i in range(2, int(D ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, D + 1, i):
                if spf[j] == j:
                    spf[j] = i
    gidx = {}                                   # d -> gpf index (1-based), squarefree p_KMAX-smooth only
    for d in range(1, D + 1):
        x, ok, g = d, True, 0
        while x > 1:
            q = spf[x]; x //= q
            if x % q == 0 or q > primes[-1]:
                ok = False; break
            g = max(g, pidx[q])
        if ok:
            gidx[d] = g
    cnt_by_g = {}
    for d, g in gidx.items():
        cnt_by_g[g] = cnt_by_g.get(g, 0) + 1
    tail_by_g = {}
    for d, g in gidx.items():
        tail_by_g[g] = tail_by_g.get(g, Fraction(0)) + Fraction(1, d ** S_TAIL)
    assert len(gidx) == J["num_small_divisors_total"], (route, len(gidx), J["num_small_divisors_total"])
    certs = {c["k"]: c for c in J["certificates"]}
    P = {s: Fraction(1) for s in J["parameters"]["s_list"]}
    Qt = Fraction(1)
    running, tsum = cnt_by_g.get(0, 0), tail_by_g.get(0, Fraction(0))
    bad = 0; worst_margin = None
    for k in range(1, KMAX + 1):
        p = primes[k - 1]
        for s in P:
            P[s] *= 1 - Fraction(1, p ** s)
        Qt *= 1 + Fraction(1, p ** S_TAIL)
        running += cnt_by_g.get(k, 0); tsum += tail_by_g.get(k, Fraction(0))
        if k < 2:
            continue
        c = certs[k]
        ok = c["num_small_divisors_used"] == running
        ep = c.get("euler_products", c.get("euler_products_exact"))
        for s in P:
            lo, hi = Fraction(ep[str(s)][0]), Fraction(ep[str(s)][1])
            ok &= lo <= P[s] <= hi
        T_exact = C_TAIL * (Qt - tsum)
        T_json = Fraction(c["tail_bound_T"])
        ok &= T_json >= T_exact                                   # their stored bound must dominate the exact tail bound
        ok &= T_json <= T_exact * (1 + Fraction(1, 10 ** 30)) + Fraction(1, 10 ** 70)   # and not be absurdly loose
        if route == "b":
            ok &= abs(Fraction(c["tail_divisor_sum_exact"][0]) - (Qt - tsum)) < Fraction(1, 10 ** 78)
        # (5) containment of my enclosure, (4) bookkeeping
        glo, ghi = Fraction(c["G_enclosure"][0]), Fraction(c["G_enclosure"][1])
        ok &= glo <= Fraction(mine[k]["G_lo"]) and Fraction(mine[k]["G_hi"]) <= ghi
        sg = -1 if ghi < 0 else (1 if glo > 0 else 0)
        ok &= sg != 0 and c["sign_G"] == sg and c["sign_a"] == (-1) ** k * sg
        pred = (k % 2 == 1 and k >= 3)
        ok &= c["conjecture_predicts_a_negative"] == pred
        ok &= c["status"].startswith("conjecture holds") == ((c["sign_a"] < 0) == pred)
        margin = min(abs(glo), abs(ghi)) / (ghi - glo)
        worst_margin = margin if worst_margin is None else min(worst_margin, margin)
        if not ok:
            bad += 1
            print("   route %s k=%d: PROBLEM  count_ok=%s euler_ok=%s tail_ok=%s contain_ok=%s sign_ok=%s" % (route, k,
                  c["num_small_divisors_used"] == running,
                  all(Fraction(ep[str(s)][0]) <= P[s] <= Fraction(ep[str(s)][1]) for s in P),
                  T_json >= T_exact, glo <= Fraction(mine[k]["G_lo"]) and Fraction(mine[k]["G_hi"]) <= ghi, sg != 0 and c["sign_G"] == sg))
        if k in (91, 92, 93):
            print("   route %s k=%d: G in [%s, %s] width %.1e, |G|/width = %.1e, sign_G=%+d sign_a=%+d status=%r, T=%.2e (exact %.2e), divisors %d"
                  % (route, k, c["G_enclosure"][0][:24], c["G_enclosure"][1][:24], float(ghi - glo), float(margin), c["sign_G"], c["sign_a"], c["status"],
                     float(T_json), float(T_exact), running))
    print("route %s (m=%d, D=%d): %d problems over k=2..%d; divisor counts, Euler products, tail bounds, bookkeeping, containment of my enclosures all verified; smallest |G|/width = %.1e"
          % (route, m, D, bad, KMAX, float(worst_margin)))
    assert bad == 0
print("cross-examination of both JSON certificates passed")
