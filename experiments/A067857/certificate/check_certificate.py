#!/usr/bin/env python3
"""
Standalone checker for the A067857 primorial sign certificates (cert_route_a.json, cert_route_b.json).

Pure Python (integers / fractions), no third-party packages.  From the JSON *inputs*
(m, D, KMAX, the prime list, the Bernoulli numbers) it RE-DERIVES a rigorous enclosure of

      G(k) = (-1)^k F(p_k#) = sum_{d | p_k#} (-1)^{omega(d)} H_d

for every k listed, using a third kind of arithmetic: fixed-point integer intervals
(scale 2^PBITS; every rounding is outward) plus exact rationals, with its own enclosures of
ln p (atanh series + tail bound) and of gamma (from H_D and the remainder lemma).  It then
verifies, for every k:
   * its own enclosure is sign-definite and has the sign claimed in the JSON;
   * the JSON enclosure and the re-derived enclosure intersect (consistency);
   * the JSON's status / sign_a fields are consistent with the sign and the parity of k;
   * the JSON's Euler-product enclosures contain the exact rational products;
   * the JSON's ln p and gamma enclosures intersect the checker's own enclosures.
Exit status 0 iff everything passes.

Mathematical basis (see README.md):
  L1: eps(d) = H_d - ln d - gamma = A_m(d) + R_m(d), A_m(d) = 1/(2d) - sum_{j<=m} B_{2j}/(2j d^{2j}),
      |R_m(d)| < |B_{2m+2}| / ((2m+2) d^{2m+2})       (DLMF 5.11(ii); Alzer 1997 Thm 8; d > 0 real)
  L2: for squarefree n with omega(n) >= 2, G = sum_{d|n} (-1)^{omega(d)} eps(d)
  L3: sum_{d|n} (-1)^{omega(d)} d^{-s} = prod_{p|n}(1-p^{-s}),  sum_{d|n} d^{-s} = prod_{p|n}(1+p^{-s})
  =>  G(k) = 1/2 P_1 - sum_j B_{2j}/(2j) P_{2j} + sum_{d<=D} (-1)^{omega(d)} (eps(d) - A_m(d)) + E,
      |E| <= |B_{2m+2}|/(2m+2) * (prod_{i<=k}(1+p_i^{-(2m+2)}) - sum_{d<=D, d|n} d^{-(2m+2)}).
Usage:  python3 check_certificate.py cert_route_a.json [cert_route_b.json ...] [--kmax K]
"""
import json, sys, time
from fractions import Fraction
from math import comb

sys.set_int_max_str_digits(0)
PBITS = 400
ONE = 1 << PBITS
JSER = 100


# ----------------------------------------------------------------- fixed-point interval helpers
def fp_from_frac(fr):
    n, d = fr.numerator, fr.denominator
    lo = (n << PBITS) // d
    hi = -((-(n << PBITS)) // d)
    return (lo, hi)

def fp_add(a, b):
    return (a[0] + b[0], a[1] + b[1])

def fp_sub(a, b):
    return (a[0] - b[1], a[1] - b[0])

def fp_neg(a):
    return (-a[1], -a[0])

def fp_scale_int(a, c):
    return (c * a[0], c * a[1]) if c >= 0 else (c * a[1], c * a[0])

def fp_to_frac(a):
    return (Fraction(a[0], ONE), Fraction(a[1], ONE))

def fp_from_frac_interval(lo, hi):
    return (fp_from_frac(lo)[0], fp_from_frac(hi)[1])


# ----------------------------------------------------------------- number theory helpers
def bernoulli_even(nmax):
    Bs = [Fraction(1)]
    for n in range(1, nmax + 1):
        Bs.append(-sum(comb(n + 1, j) * Bs[j] for j in range(n)) / (n + 1))
    return {2 * j: Bs[2 * j] for j in range(1, nmax // 2 + 1)}

def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True

def atanh_enclosure(t):
    assert 0 <= t <= Fraction(1, 3)
    t2 = t * t
    term = t
    s = Fraction(0)
    for j in range(JSER):
        s += term / (2 * j + 1)
        term *= t2
    return s, s + term / ((2 * JSER + 1) * (1 - t2))

LN2 = tuple(2 * x for x in atanh_enclosure(Fraction(1, 3)))

def ln_int_enclosure(p):
    if p == 1:
        return (Fraction(0), Fraction(0))
    e = p.bit_length() - 1
    lo, hi = atanh_enclosure(Fraction(p - 2 ** e, p + 2 ** e))
    return (e * LN2[0] + 2 * lo, e * LN2[1] + 2 * hi)

def parse_enc(pair):
    lo, hi = Fraction(pair[0]), Fraction(pair[1])
    assert lo <= hi
    return lo, hi

def overlaps(a, b):
    return a[0] <= b[1] and b[0] <= a[1]


# ----------------------------------------------------------------- main verification of one JSON
def verify(path, kmax_limit=None):
    t0 = time.time()
    J = json.load(open(path))
    par = J["parameters"]
    m, D, KMAX = par["m"], par["D"], par["KMAX"]
    primes = J["primes"]
    if kmax_limit is not None:
        KMAX = min(KMAX, kmax_limit)
    primes = primes[:KMAX]
    fails = []
    def check(cond, msg):
        if not cond:
            fails.append(msg)
            print("   FAIL:", msg)

    print("== %s : route %s ; m=%d D=%d KMAX=%d" % (path, J.get("route", "?"), m, D, KMAX))
    # -- primes: the first KMAX primes, verified independently
    q, cnt, mine = 2, 0, []
    while cnt < KMAX:
        if is_prime(q):
            mine.append(q); cnt += 1
        q += 1
    check(mine == primes, "prime list is not the first %d primes" % KMAX)

    # -- Bernoulli numbers and coefficients
    B = bernoulli_even(2 * m + 2)
    for s, v in J["bernoulli"].items():
        check(Fraction(v) == B[int(s)], "Bernoulli number B_%s mismatch" % s)
    coef = {1: Fraction(1, 2)}
    for j in range(1, m + 1):
        coef[2 * j] = -B[2 * j] / (2 * j)
    for s, v in J["coefficients"].items():
        check(Fraction(v) == coef[int(s)], "coefficient for s=%s mismatch" % s)
    S_TAIL = 2 * m + 2
    C_TAIL = abs(B[S_TAIL]) / S_TAIL
    check(Fraction(J["tail_constant"]) == C_TAIL, "tail constant mismatch")
    S_LIST = sorted(coef)

    def A_m(d):
        d = Fraction(d)
        return 1 / (2 * d) - sum(B[2 * j] / (2 * j * d ** (2 * j)) for j in range(1, m + 1))

    # -- own ln p enclosures; compare with the JSON's
    LN = {p: ln_int_enclosure(p) for p in primes}
    for p in primes:
        check(overlaps(LN[p], parse_enc(J["constants"]["ln_p"][str(p)])), "ln %d enclosure inconsistent" % p)
    LNfp = {p: fp_from_frac_interval(*LN[p]) for p in primes}

    # -- divisors d <= D of p_KMAX# (squarefree, primes among the first KMAX)
    divs = []
    def dfs(start, d, idxs):
        divs.append((d, len(idxs), (idxs[-1] + 1) if idxs else 0, tuple(idxs)))
        for i in range(start, KMAX):
            nd = d * primes[i]
            if nd > D:
                break
            dfs(i + 1, nd, idxs + [i])
    dfs(0, 1, [])
    div_info = {d: (w, j, idxs) for d, w, j, idxs in divs}
    # sanity: every squarefree p_KMAX-smooth integer <= D is in the list (brute force check for D <= 10^5)
    if D <= 10 ** 5:
        cntchk = 0
        pset = set(primes)
        for x in range(1, D + 1):
            y, ok = x, True
            for p in primes:
                if p * p > y:
                    break
                if y % p == 0:
                    y //= p
                    if y % p == 0:
                        ok = False; break
            if ok and (y == 1 or y in pset):
                cntchk += 1
        check(cntchk == len(divs), "divisor enumeration count mismatch (%d vs %d)" % (cntchk, len(divs)))

    # -- fixed-point interval harmonic numbers H_d, d <= D
    Hfp = {}
    lo = hi = 0
    for i in range(1, D + 1):
        lo += ONE // i
        hi += -((-ONE) // i)
        if i in div_info or i == D:
            Hfp[i] = (lo, hi)
    HDfp = Hfp[D]

    # -- gamma from H_D: gamma = H_D - ln D - A_m(D) - R_m(D), |R_m(D)| < C_TAIL / D^{2m+2}
    # ln D: D must be 10^a here
    a10 = len(str(D)) - 1
    check(D == 10 ** a10, "checker expects D to be a power of 10")
    lnD = fp_from_frac_interval(a10 * (LN[2][0] + LN[5][0]), a10 * (LN[2][1] + LN[5][1]))
    RD = C_TAIL / Fraction(D) ** S_TAIL
    gam = fp_sub(fp_sub(fp_sub(HDfp, lnD), fp_from_frac(A_m(D))), fp_from_frac_interval(-RD, RD))
    gam_fr = fp_to_frac(gam)
    check(overlaps(gam_fr, parse_enc(J["constants"]["gamma"])), "gamma enclosure inconsistent with JSON")
    check(Fraction("0.57721566490153286060651209008240243104215") < gam_fr[0]
          and gam_fr[1] < Fraction("0.57721566490153286060651209008240243104216"), "gamma enclosure off (vs known digits)")
    print("   own gamma enclosure width %.2e ; ln p widths <= %.2e ; %d small divisors ; H_D width %.2e   [%.1fs]"
          % (float(gam_fr[1] - gam_fr[0]), max(float(LN[p][1] - LN[p][0]) for p in primes), len(divs),
             float(Fraction(HDfp[1] - HDfp[0], ONE)), time.time() - t0))

    # -- per-divisor corrections v(d) = (-1)^omega (H_d - ln d - gamma - A_m(d)), bucketed by gpf index
    bucketV = [(0, 0)] * (KMAX + 1)
    bucketS = [(0, 0)] * (KMAX + 1)
    for d, w, j, idxs in divs:
        v = Hfp[d]
        for i in idxs:
            v = fp_sub(v, LNfp[primes[i]])
        v = fp_sub(v, gam)
        v = fp_sub(v, fp_from_frac(A_m(d)))
        if w % 2:
            v = fp_neg(v)
        bucketV[j] = fp_add(bucketV[j], v)
        bucketS[j] = fp_add(bucketS[j], fp_from_frac(Fraction(1, d ** S_TAIL)))

    # -- per k
    P = {s: Fraction(1) for s in S_LIST}
    Qt = Fraction(1)
    corr = bucketV[0]
    ssum = bucketS[0]
    certs = {c["k"]: c for c in J["certificates"]}
    n_ok = 0
    first_neg = None
    signs = {}
    for k in range(1, KMAX + 1):
        p = primes[k - 1]
        for s in S_LIST:
            P[s] *= 1 - Fraction(1, p ** s)
        Qt *= 1 + Fraction(1, p ** S_TAIL)
        corr = fp_add(corr, bucketV[k])
        ssum = fp_add(ssum, bucketS[k])
        if k < 2:
            continue
        closed = sum(coef[s] * P[s] for s in S_LIST)
        tail_rest_hi = max(Fraction(0), Qt - Fraction(ssum[0], ONE))       # >= sum_{d|n, d>D} d^{-(2m+2)}
        T = C_TAIL * tail_rest_hi
        G = fp_add(fp_add(fp_from_frac(closed), corr), fp_from_frac_interval(-T, T))
        Glo, Ghi = fp_to_frac(G)
        sgn = -1 if Ghi < 0 else (1 if Glo > 0 else 0)
        signs[k] = sgn
        if sgn < 0 and first_neg is None:
            first_neg = k
        c = certs.get(k)
        if c is None:
            check(False, "k=%d missing from JSON" % k)
            continue
        jlo, jhi = parse_enc(c["G_enclosure"])
        ok = True
        ok &= sgn != 0
        ok &= (c["sign_G"] == sgn)
        ok &= overlaps((Glo, Ghi), (jlo, jhi))
        ok &= (c["sign_a"] == (-1) ** k * sgn)
        pred = (k % 2 == 1 and k >= 3)
        ok &= (c["conjecture_predicts_a_negative"] == pred)
        holds = ((c["sign_a"] < 0) == pred)
        ok &= (c["status"].startswith("conjecture holds") == holds)
        # Euler products in the JSON must contain the exact rational values
        ep = c.get("euler_products", c.get("euler_products_exact"))
        for s in S_LIST:
            elo, ehi = parse_enc(ep[str(s)])
            ok &= (elo <= P[s] <= ehi)
        if not ok:
            check(False, "k=%d: re-derived G in [%s, %s] (sign %+d) vs JSON [%s, %s] sign %+d status %r"
                  % (k, float(Glo), float(Ghi), sgn, c["G_enclosure"][0][:25], c["G_enclosure"][1][:25], c["sign_G"], c["status"]))
        else:
            n_ok += 1
        if k in (2, 3, 88, 89, 90, 91, 92, 93, 94, 95, 96, KMAX):
            print("   k=%3d p_k=%4d  re-derived G in [%+.9e, %+.9e] (width %.1e)  sign %+d  JSON sign %+d  %s"
                  % (k, p, float(Glo), float(Ghi), float(Ghi - Glo), sgn, c["sign_G"], "ok" if ok else "MISMATCH"))
    check(J["first_k_with_G_negative"] == first_neg, "first negative k mismatch (JSON %s, checker %s)" % (J["first_k_with_G_negative"], first_neg))
    pos = [k for k in signs if signs[k] > 0]
    neg = [k for k in signs if signs[k] < 0]
    und = [k for k in signs if signs[k] == 0]
    print("   re-derived: G(k) > 0 for k in [%s..%s] (%d values), G(k) < 0 for k in [%s..%s] (%d values), undecided: %s"
          % (min(pos) if pos else None, max(pos) if pos else None, len(pos), min(neg) if neg else None, max(neg) if neg else None, len(neg), und))
    check(pos == list(range(2, first_neg)) and neg == list(range(first_neg, KMAX + 1)) and not und,
          "sign pattern is not 'positive then negative'")
    print("   %d/%d certificates verified   [%.1fs]" % (n_ok, KMAX - 1, time.time() - t0))
    return fails


if __name__ == "__main__":
    argv = sys.argv[1:]
    kmax_limit = None
    if "--kmax" in argv:
        i = argv.index("--kmax")
        kmax_limit = int(argv[i + 1])
        del argv[i:i + 2]
    args = argv
    if not args:
        args = ["cert_route_a.json", "cert_route_b.json"]
    allfails = []
    for path in args:
        allfails += verify(path, kmax_limit)
    if allfails:
        print("CHECK FAILED: %d problem(s)" % len(allfails))
        sys.exit(1)
    print("ALL CHECKS PASSED")
