"""lemmas-audit: independent exact counter + tests of identity (i) and of the derived SG caps (ii)
on real point sets.  Written from scratch (no import of the theory agent's code).

Block key: a circle is keyed by (2*cx, 2*cy, r^2) computed with Fractions; a line by its primitive
integer equation.  A block = circle-or-line through >= 3 points.
"""
from fractions import Fraction as Fr
from itertools import combinations
from math import comb, gcd
import random, sys

def f(n):
    return comb(n - 1, 2) + 1 - (n - 1) // 2

def line_key(p, q):
    (x1, y1), (x2, y2) = p, q
    a, b = y2 - y1, x1 - x2
    c = a * x1 + b * y1                     # a x + b y = c
    den = 1
    for v in (a, b, c):
        den = den * v.denominator // gcd(den, v.denominator)
    a, b, c = int(a * den), int(b * den), int(c * den)
    g = gcd(gcd(abs(a), abs(b)), abs(c))
    a, b, c = a // g, b // g, c // g
    if a < 0 or (a == 0 and b < 0):
        a, b, c = -a, -b, -c
    return ("L", a, b, c)

def circle_key(p, q, r):
    (x1, y1), (x2, y2), (x3, y3) = p, q, r
    det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if det == 0:
        return None
    s1, s2, s3 = x1 * x1 + y1 * y1, x2 * x2 + y2 * y2, x3 * x3 + y3 * y3
    # solve 2(x2-x1)u + 2(y2-y1)v = s2-s1 ; 2(x3-x1)u + 2(y3-y1)v = s3-s1
    u = ((s2 - s1) * (y3 - y1) - (s3 - s1) * (y2 - y1)) / (2 * det)
    v = ((x2 - x1) * (s3 - s1) - (x3 - x1) * (s2 - s1)) / (2 * det)
    return ("C", u, v, (x1 - u) ** 2 + (y1 - v) ** 2)

def blocks(points):
    pts = [(Fr(x), Fr(y)) for x, y in points]
    assert len(set(pts)) == len(pts)
    B = {}
    for i, j, k in combinations(range(len(pts)), 3):
        key = circle_key(pts[i], pts[j], pts[k])
        if key is None:
            key = line_key(pts[i], pts[j])
        B.setdefault(key, set()).update((i, j, k))
    return B

def analyse(points):
    n = len(points)
    B = blocks(points)
    circles = [frozenset(s) for k, s in B.items() if k[0] == "C"]
    lines = [frozenset(s) for k, s in B.items() if k[0] == "L"]
    # sanity: every triple in exactly one block, blocks pairwise share <= 2 points, lines <= 1
    allb = circles + lines
    assert sum(comb(len(b), 3) for b in allb) == comb(n, 3)
    assert all(len(a & b) <= 2 for a, b in combinations(allb, 2))
    assert all(len(a & b) <= 1 for a, b in combinations(lines, 2))
    D = sum(comb(len(b), 3) - 1 for b in allb if len(b) >= 4)
    return dict(n=n, circles=len(circles), lines=len(lines), D=D,
                identity=(len(circles) == comb(n, 3) - D - len(lines)),
                degenerate=any(len(b) == n for b in allb),
                csizes=sorted(len(b) for b in circles), lsizes=sorted(len(b) for b in lines),
                blocks=allb, line_blocks=lines)

def derived_cap_check(points, res):
    """(ii): for each p, pairs of P\\{p} covered by blocks of size>=4 through p <= C(n-1,2)-1
    (SG), and more precisely the number of ordinary derived lines = pairs {a,b} whose block with p
    has size 3, must be >= 1; lines: C(n,2) - sum C(|l|,2) >= 1."""
    n = res["n"]
    out = []
    for p in range(n):
        thru = [b for b in res["blocks"] if p in b and len(b) >= 4]
        covered = sum(comb(len(b) - 1, 2) for b in thru)
        ordinary = comb(n - 1, 2) - covered
        out.append(ordinary)
    lines_ord = comb(n, 2) - sum(comb(len(l), 2) for l in res["line_blocks"])
    return out, lines_ord

RECORDS = {
    "n=6 (report)": [(0, 0), (20, 0), (5, 15), (5, 0), (2, 6), (10, 10)],
    "n=7 (report)": [(0, 0), (20, 0), (5, 15), (5, 0), (2, 6), (10, 10), (5, 5)],
    "n=8 (report)": [(0, 0), (0, 5), (0, 10), (0, 15), (3, 6), (3, 9), (5, 5), (5, 10)],
    "n=8 two squares": [(1, 0), (0, 1), (-1, 0), (0, -1), (2, 0), (0, 2), (-2, 0), (0, -2)],
    "n=9 antipodal": [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1), (Fr(3, 5), Fr(4, 5)),
                       (Fr(-3, 5), Fr(-4, 5)), (Fr(3, 5), Fr(-4, 5)), (Fr(-3, 5), Fr(4, 5))],
}

def rational_circle_points(cx, cy, r, ts):
    """rational points on circle centre (cx,cy) radius r via t -> ((1-t^2)/(1+t^2), 2t/(1+t^2))."""
    out = []
    for t in ts:
        t = Fr(t)
        out.append((cx + r * (1 - t * t) / (1 + t * t), cy + r * 2 * t / (1 + t * t)))
    return out

def random_structured_set(rng, n):
    """random rational set with forced concyclic / collinear subsets, to test identity (i)."""
    pts = set()
    while len(pts) < n:
        kind = rng.random()
        if kind < 0.4:      # points on a random circle
            cx, cy, r = Fr(rng.randint(-3, 3)), Fr(rng.randint(-3, 3)), Fr(rng.randint(1, 3))
            for p in rational_circle_points(cx, cy, r, [Fr(rng.randint(-5, 5), rng.randint(1, 4)) for _ in range(rng.randint(3, 5))]):
                pts.add(p)
        elif kind < 0.7:    # points on a random line
            x0, y0 = Fr(rng.randint(-3, 3)), Fr(rng.randint(-3, 3))
            dx, dy = Fr(rng.randint(-3, 3)), Fr(rng.randint(-3, 3))
            if dx == 0 and dy == 0:
                continue
            for _ in range(rng.randint(3, 4)):
                t = Fr(rng.randint(-4, 4), rng.randint(1, 3))
                pts.add((x0 + t * dx, y0 + t * dy))
        else:
            pts.add((Fr(rng.randint(-6, 6), rng.randint(1, 3)), Fr(rng.randint(-6, 6), rng.randint(1, 3))))
    pts = list(pts)[:n]
    return pts

if __name__ == "__main__":
    print("== constructions ==")
    for name, P in RECORDS.items():
        r = analyse(P)
        ords, lo = derived_cap_check(P, r)
        print(f"{name}: n={r['n']} circles={r['circles']} f(n)={f(r['n'])} lines={r['lines']} D={r['D']} "
              f"identity={r['identity']} degenerate={r['degenerate']} csizes={r['csizes']} lsizes={r['lsizes']}")
        print(f"    ordinary derived lines per point (must all be >= 1): {ords}; ordinary lines at infinity: {lo}")
    print("== identity (i) on random structured sets ==")
    rng = random.Random(12345)
    bad = 0; tested = 0
    for trial in range(400):
        n = rng.randint(4, 10)
        P = random_structured_set(rng, n)
        if len(P) < 4:
            continue
        r = analyse(P)
        if r["degenerate"]:
            continue
        tested += 1
        ords, lo = derived_cap_check(P, r)
        if not r["identity"] or min(ords) < 1 or lo < 1:
            bad += 1
            print("VIOLATION", P, r)
    print(f"tested {tested} non-degenerate random sets; identity/cap violations: {bad}")
