#!/usr/bin/env python3
"""postcheck.py -- hereditary realisability checks for candidate structures produced by orderly2.py.

For each candidate (rich blocks F on n points, list of line sets L achieving the required number of lines):
 * derived structure at p: the family {B \\ {p} : p in B in F} of >=3-point lines of P \\ {p}
   (after inversion about p these are exactly the rich lines of a real non-collinear (n-1)-point set);
 * line structure at infinity: L itself is the family of >=3-point lines of the real n-point set P.
Tests applied to every such linear space (X, lines):
 (H1) hereditary Sylvester-Gallai: for every S subset X with |S| >= 3 not contained in one line, some pair
      of S lies in no line containing a third point of S.  A violating S is "SG-closed"; the Fano plane and
      AG(2,3) are the classical examples.
 (H2) hereditary orchard: for every S, #{lines meeting S in exactly 3 points} <= t3(|S|),
      t3(3..12) = 1,1,2,4,6,7,10,12,16,19 (Burr-Gruenbaum-Sloane 1974, OEIS A003035).  In particular a
      Moebius-Kantor (8_3) sub-configuration (8 three-point lines on 8 points) violates t3(8) = 7.
 (H3) ordinary-line count: the number of pairs of X lying in no line is >= o(|X|) (only for the full X;
      this is the cap already used in the enumeration, repeated here for the record).
A candidate survives iff every derived structure passes (H1)-(H3) and at least one of its line sets passes
(H1)-(H3) at infinity.
Usage: python3 postcheck.py results.json
"""
import sys, json, itertools
from math import comb

T3 = {3: 1, 4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10, 10: 12, 11: 16, 12: 19}
O_TABLE = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5, 11: 6, 12: 6, 13: 6, 14: 7}


def check_linear_space(X, lines, use_table=True):
    """X: list of points; lines: list of frozensets (each of size >= 3, pairwise sharing <= 1 point).
    Returns a list of violation descriptions (empty = passes)."""
    X = sorted(X)
    viol = []
    # (H3)
    covered = set()
    for l in lines:
        for pr in itertools.combinations(sorted(l), 2):
            covered.add(pr)
    ordinary = comb(len(X), 2) - len(covered)
    need = O_TABLE[len(X)] if use_table else 1
    if len(X) >= 3 and not any(set(X) <= l for l in lines) and ordinary < need:
        viol.append(f"H3: only {ordinary} ordinary lines on {len(X)} points (need >= {need})")
    # (H1), (H2) over all subsets
    for r in range(3, len(X) + 1):
        for S in itertools.combinations(X, r):
            Sset = set(S)
            if any(Sset <= l for l in lines):
                continue
            restricted = [l & Sset for l in lines if len(l & Sset) >= 3]
            cov = set()
            for l in restricted:
                for pr in itertools.combinations(sorted(l), 2):
                    cov.add(pr)
            if len(cov) == r * (r - 1) // 2:
                viol.append(f"H1: SG-closed subset {list(S)} with lines {[sorted(l) for l in restricted]}")
            n3 = sum(1 for l in restricted if len(l) == 3)
            if n3 > T3[r]:
                viol.append(f"H2: subset {list(S)} has {n3} three-point lines > t3({r}) = {T3[r]}")
    return viol


def check_candidate(rec, use_table=True):
    n = rec["n"]
    blocks = [frozenset(b) for b in rec["blocks"]]
    report = {"derived": {}, "lines": []}
    for p in range(n):
        derived = [b - {p} for b in blocks if p in b]
        v = check_linear_space([q for q in range(n) if q != p], derived, use_table)
        if v:
            report["derived"][p] = v
    for L in rec.get("line_sets", []):
        v = check_linear_space(list(range(n)), [frozenset(l) for l in L], use_table)
        report["lines"].append(v)
    survives = (not report["derived"]) and any(not v for v in report["lines"]) if rec.get("line_sets") else \
        (not report["derived"])
    return survives, report


if __name__ == "__main__":
    data = json.load(open(sys.argv[1]))
    recs = data["results"] if isinstance(data, dict) else data
    print(f"{len(recs)} candidate structures")
    surv = 0
    for i, rec in enumerate(recs):
        ok, rep = check_candidate(rec)
        print(f"structure {i}: sizes={rec['sizes']} D={rec['D']} ell_max={rec['ell_max']} count_min={rec['count_min']} "
              f"degrees={rec['degrees']} -> {'SURVIVES' if ok else 'EXCLUDED'}")
        print("   blocks:", rec["blocks"])
        for p, v in rep["derived"].items():
            print(f"   derived structure at point {p}: {len(v)} violations, e.g. {v[0]}")
        bad = [j for j, v in enumerate(rep["lines"]) if v]
        print(f"   line sets: {len(rep['lines'])} total, {len(bad)} violate; "
              f"{'e.g. ' + rep['lines'][bad[0]][0] if bad else ''}")
        surv += ok
    print(f"survivors: {surv}")
