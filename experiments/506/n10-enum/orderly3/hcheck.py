#!/usr/bin/env python3
"""
hcheck.py -- hereditary realisability checks for candidate structures (F, L) produced by oe.py.

For the 11-point Moebius set Q = P + {inf} every point q of Q has a *derived* real point set
Q - {q} (obtained by a Moebius transformation sending q to infinity) whose lines with >= 3 points
are exactly the blocks of Q through q, minus q.  Concretely, for p in P:
    derived lines at p  = { B - p             : B in F, p in B, B not a line }
                        + { (B - p) + {inf'}  : B in F, p in B, B in L }
                        + { (T - p) + {inf'}  : T a 3-line of L, p in T }
    (points: P - p and inf' = image of inf, 10 points)
and for q = inf the derived structure is the line set L itself on P.
The 9-point structure "F at p" (blocks through p, no inf) is the restriction to P - p.

Checks on a real point set X with known complete list R of its >= 3-point lines:
  tier A (theorem: Sylvester-Gallai, hereditary): for every S subset of X, |S| >= 3, S not inside one
         line: some pair of S is on no line of R with a third point of S.
         (A 7-subset carrying a Fano plane is the classical violation.)
  tier B (theorem proved in this session, ../../verify_independent/mk_unrealisable2.py): no 8-subset S
         whose >= 3-point lines (restricted) are exactly 8 triples pairwise sharing <= 1 point
         (Moebius-Kantor 8_3).
  tier C (cited tables): for every S: #(restricted lines with exactly 3 points of S) <= t3(|S|)
         (Burr-Gruenbaum-Sloane orchard numbers) and #ordinary lines of S >= o(|S|) if S is not
         collinear (Kelly-Moser / Csima-Sawyer / Crowe-McKee values).
A structure F survives a tier if at least one of its line sets L (with count <= target) passes all
checks of that tier at all 11 points of Q.
"""
import sys, json, itertools
from math import comb
from oe import O_TABLE, T3_TABLE, popcount, mask_to_list, list_to_mask


def restricted_lines(R, S):
    """R: list of line masks; S: subset mask.  Lines of R with >= 3 points in S, restricted."""
    out = []
    for L in R:
        x = L & S
        if popcount(x) >= 3:
            out.append(x)
    return out


def check_pointset(m, R, tiers=('A', 'B', 'K', 'C')):
    """m points 0..m-1, R = complete list of >= 3-point lines (masks).  Returns dict tier -> list of
    violation strings (empty list = passes)."""
    viol = {t: [] for t in tiers}
    full = (1 << m) - 1
    for S in range(1, 1 << m):
        s = popcount(S)
        if s < 3:
            continue
        RS = restricted_lines(R, S)
        if any(x == S for x in RS):          # S collinear
            continue
        covered = sum(comb(popcount(x), 2) for x in RS)
        ordinary = comb(s, 2) - covered
        if 'A' in tiers and ordinary == 0:
            viol['A'].append(f"SG-closed subset {mask_to_list(S, m)} lines {[mask_to_list(x, m) for x in RS]}")
        if 'B' in tiers and s == 8:
            n3 = sum(1 for x in RS if popcount(x) == 3)
            if n3 == len(RS) and n3 >= 8:
                viol['B'].append(f"(8_3) on {mask_to_list(S, m)} lines {[mask_to_list(x, m) for x in RS]}")
        if 'K' in tiers and s >= 3 and ordinary < -(-3 * s // 7):     # Kelly-Moser: o(m) >= 3m/7
            viol['K'].append(f"subset {mask_to_list(S, m)} has {ordinary} ordinary lines < ceil(3*{s}/7)={-(-3 * s // 7)}")
        if 'C' in tiers:
            n3 = sum(1 for x in RS if popcount(x) == 3)
            if s in T3_TABLE and n3 > T3_TABLE[s]:
                viol['C'].append(f"subset {mask_to_list(S, m)} has {n3} three-point lines > t3({s})={T3_TABLE[s]}")
            if s in O_TABLE and ordinary < O_TABLE[s]:
                viol['C'].append(f"subset {mask_to_list(S, m)} has {ordinary} ordinary lines < o({s})={O_TABLE[s]}")
    return viol


def derived_at_p(n, F, L, p):
    """Derived 10-point structure at p of Q = P + inf (inf' gets index n).  Returns line masks."""
    Lset = set(L)
    R = []
    for B in F:
        if B >> p & 1:
            x = B & ~(1 << p)
            if B in Lset:
                x |= 1 << n
            R.append(x)
    for T in L:
        if popcount(T) == 3 and T >> p & 1:
            R.append((T & ~(1 << p)) | (1 << n))
    return R


def relabel(mask, pts):
    """restrict mask to the point list pts and relabel to 0..len(pts)-1."""
    return sum(1 << j for j, q in enumerate(pts) if mask >> q & 1)


def check_structure(n, blocks, line_sets, tiers=('A', 'B', 'K', 'C'), verbose=False):
    F = [list_to_mask(b) for b in blocks]
    report = {"F_only": {}, "line_sets": []}
    # 9-point derived structures of F alone (independent of L)
    f_viol = {t: [] for t in tiers}
    for p in range(n):
        pts = [q for q in range(n) if q != p]
        R = [relabel(B, pts) for B in F if B >> p & 1]
        v = check_pointset(n - 1, R, tiers)
        for t in tiers:
            if v[t]:
                f_viol[t].append((p, v[t][0], len(v[t])))
    report["F_only"] = f_viol
    survivors = {t: [] for t in tiers}
    for li, Ls in enumerate(line_sets):
        L = [list_to_mask(l) for l in Ls]
        lv = {t: [] for t in tiers}
        # at infinity: the line set itself on P
        v = check_pointset(n, L, tiers)
        for t in tiers:
            if v[t]:
                lv[t].append(("inf", v[t][0], len(v[t])))
        # at each p: 10-point derived structure with inf'
        for p in range(n):
            pts = [q for q in range(n) if q != p] + [n]
            R = [relabel(x, pts) for x in derived_at_p(n, F, L, p)]
            v = check_pointset(n, R, tiers)
            for t in tiers:
                if v[t]:
                    lv[t].append((p, v[t][0], len(v[t])))
        report["line_sets"].append(lv)
        for t in tiers:
            if not lv[t] and not f_viol[t]:
                survivors[t].append(li)
    report["survivors"] = survivors
    return report


def cumulative(report, tiers=('A', 'B', 'K', 'C')):
    """line sets surviving tier A, tiers A+B, tiers A+B+C."""
    s = report["survivors"]
    out = {}
    acc = None
    for t in tiers:
        acc = set(s[t]) if acc is None else acc & set(s[t])
        out[t] = sorted(acc)
    return out


def main():
    fn = sys.argv[1]
    data = json.load(open(fn))
    n = data["n"]
    results = data["results"]
    print(f"{fn}: {len(results)} candidate structures (n={n}, target={data['target']}, mode={data['mode']})")
    summary = []
    for i, r in enumerate(results):
        rep = check_structure(n, r["blocks"], r["line_sets"])
        cum = cumulative(rep)
        counts = {t: sorted(set(comb(n, 3) - r["D"] - len(r["line_sets"][j]) for j in cum[t])) for t in cum}
        print(f"structure {i}: sizes={r['sizes']} D={r['D']} ell_max={r['ell_max']} count_min={r['count_min']} "
              f"degrees={r['degrees']} #line_sets={len(r['line_sets'])}")
        print(f"   blocks={r['blocks']}")
        for t in ('A', 'B', 'K', 'C'):
            fv = rep['F_only'][t]
            if fv:
                p, ex, k = fv[0]
                print(f"   F alone, tier {t}: violated at {len(fv)} points, e.g. p={p}: {ex}")
        print(f"   line sets surviving: tier A: {len(cum['A'])} (counts {counts['A']}); "
              f"A+B: {len(cum['B'])} (counts {counts['B']}); A+B+K: {len(cum['K'])} (counts {counts['K']}); "
              f"A+B+K+C: {len(cum['C'])} (counts {counts['C']})")
        if rep['line_sets']:
            # show one example violation for the first line set for each tier
            lv = rep['line_sets'][0]
            for t in ('A', 'B', 'K', 'C'):
                if lv[t]:
                    q, ex, k = lv[t][0]
                    print(f"   line set 0, tier {t}: e.g. at {q}: {ex}")
        verdict = ("SURVIVES all tiers" if cum['C'] else
                   "killed by tier C only (orchard / exact o-table)" if cum['K'] else
                   "killed by tier K (Kelly-Moser o(m) >= 3m/7)" if cum['B'] else
                   "killed by tier B ((8_3) lemma)" if cum['A'] else "killed by tier A (hereditary SG)")
        print(f"   => {verdict}")
        summary.append((i, r['sizes'].count(5), r['sizes'].count(6), len(r['blocks']), r['count_min'],
                        len(cum['A']), len(cum['B']), len(cum['K']), len(cum['C']), counts, verdict))
    print("\nSUMMARY (structure, #5-blocks, #6-blocks, #blocks, count_min, surviving line sets A / A+B / A+B+K / A+B+K+C, counts, verdict)")
    for s in summary:
        print("  ", s)
    if len(sys.argv) > 2:
        json.dump([{"index": s[0], "verdict": s[10], "surv_A": s[5], "surv_AB": s[6], "surv_ABK": s[7],
                    "surv_ABKC": s[8], "counts": {k: v for k, v in s[9].items()}} for s in summary],
                  open(sys.argv[2], 'w'), indent=1)


if __name__ == '__main__':
    main()
