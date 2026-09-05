"""Hereditary Sylvester-Gallai checks for candidate structures produced by mobius_enum*.py.

For a structure (blocks F on n points) and a point p, the derived linear space at p is the family
{B \ {p} : p in B, |B| >= 4} on P \ {p}; after inversion about p these are exactly the >=3-point
lines of a real (n-1)-point set.  Sylvester-Gallai applied to every subset S of that point set says:
if S is not contained in one derived line then some pair of S lies on no derived line together with
a third point of S.  Any violating subset S (a 'SG-closed' subset, e.g. a Fano plane) proves the
structure non-realisable.  The same test applies to the lines L at the point 'infinity'.
Usage: python3 hereditary_sg.py enum_output.json
"""
import sys, json, itertools

def sg_violations(points, lines):
    """lines: list of frozensets on 'points'.  Return a list of subsets S (as sorted lists) that are
    SG-closed: |S|>=3, S not inside one line, every pair of S lies in a line whose intersection with
    S has >= 3 points.  Only minimal violators (no proper violating subset) are returned."""
    pts = sorted(points)
    viol = []
    for r in range(3, len(pts) + 1):
        for S in itertools.combinations(pts, r):
            Sset = set(S)
            if any(Sset <= l for l in lines):
                continue
            restricted = [l & Sset for l in lines if len(l & Sset) >= 3]
            covered = set()
            for l in restricted:
                for pair in itertools.combinations(sorted(l), 2):
                    covered.add(pair)
            if len(covered) == r * (r - 1) // 2:
                if not any(set(v) <= Sset for v in viol):
                    viol.append(list(S))
    return viol

def check_structure(rec):
    n = rec["n"]
    blocks = [frozenset(b) for b in rec["blocks"]]
    report = {}
    for p in range(n):
        derived = [b - {p} for b in blocks if p in b]
        v = sg_violations(set(range(n)) - {p}, derived)
        if v:
            report[p] = v
    # lines at infinity: each line set
    line_report = []
    for L in rec.get("line_sets", []):
        v = sg_violations(set(range(n)), [frozenset(l) for l in L])
        line_report.append(v)
    return report, line_report

if __name__ == "__main__":
    recs = json.load(open(sys.argv[1]))
    for i, rec in enumerate(recs):
        rep, lrep = check_structure(rec)
        print(f"structure {i}: sizes={rec['sizes']} count_min={rec['count_min']} degrees={rec['degrees']}")
        if rep:
            for p, v in rep.items():
                print(f"   derived structure at point {p} violates hereditary SG on subsets: {v}")
        else:
            print("   no hereditary-SG violation in derived structures")
        bad = [j for j, v in enumerate(lrep) if v]
        print(f"   line sets violating hereditary SG: {len(bad)} of {len(lrep)}")
