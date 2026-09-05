"""Independent exact counter for Erdos #506 (audit of search-local).

A block through three integer points p1,p2,p3 is the unique curve
    A (x^2+y^2) + B x + C y + D = 0
with (A,B,C,D) obtained from the 4x4 determinant expansion; A = 0 iff collinear.
The tuple is normalised (gcd 1, first nonzero positive) so that equal blocks give equal keys.
All arithmetic is in Python integers (no Fractions, no floats).
"""
from itertools import combinations
from math import gcd, comb
from collections import defaultdict


def block_key(p, q, r):
    (x1, y1), (x2, y2), (x3, y3) = p, q, r
    s1, s2, s3 = x1 * x1 + y1 * y1, x2 * x2 + y2 * y2, x3 * x3 + y3 * y3
    # determinant | x^2+y^2  x  y  1 | expanded along first row
    A = x1 * (y2 - y3) - y1 * (x2 - x3) + (x2 * y3 - x3 * y2)
    B = -(s1 * (y2 - y3) - y1 * (s2 - s3) + (s2 * y3 - s3 * y2))
    C = s1 * (x2 - x3) - x1 * (s2 - s3) + (s2 * x3 - s3 * x2)
    D = -(s1 * (x2 * y3 - x3 * y2) - x1 * (s2 * y3 - s3 * y2) + y1 * (s2 * x3 - s3 * x2))
    g = 0
    for v in (A, B, C, D):
        g = gcd(g, abs(v))
    if g == 0:
        raise ValueError("degenerate triple (duplicate point?)")
    A, B, C, D = A // g, B // g, C // g, D // g
    for v in (A, B, C, D):
        if v != 0:
            if v < 0:
                A, B, C, D = -A, -B, -C, -D
            break
    return (A, B, C, D)


def on_block(key, pt):
    A, B, C, D = key
    x, y = pt
    return A * (x * x + y * y) + B * x + C * y + D == 0


def analyse(points):
    pts = [(int(x), int(y)) for x, y in points]
    n = len(pts)
    assert len(set(pts)) == n, "duplicate points"
    blocks = defaultdict(set)
    for i, j, k in combinations(range(n), 3):
        blocks[block_key(pts[i], pts[j], pts[k])].update((i, j, k))
    # consistency: every point on the curve is in the block, block size k has C(k,3) triples
    ntrip = 0
    for key, members in blocks.items():
        for idx, pt in enumerate(pts):
            if on_block(key, pt):
                assert idx in members, ("incidence mismatch", key, pt)
        ntrip += comb(len(members), 3)
    assert ntrip == comb(n, 3)
    circles = {k: v for k, v in blocks.items() if k[0] != 0}
    lines = {k: v for k, v in blocks.items() if k[0] == 0}
    return {
        "n": n,
        "circles": len(circles),
        "lines": len(lines),
        "circle_sizes": sorted((len(v) for v in circles.values()), reverse=True),
        "line_sizes": sorted((len(v) for v in lines.values()), reverse=True),
        "all_collinear": any(len(v) == n for v in lines.values()),
        "all_concyclic": any(len(v) == n for v in circles.values()),
        "blocks": blocks,
    }


def f(n):
    return comb(n - 1, 2) + 1 - (n - 1) // 2


if __name__ == "__main__":
    import sys, json
    from collections import Counter
    # records as given in the search-local report (task statement)
    REC = {
        6: (8, "(15,0),(0,15),(10,5),(0,-5),(6,-3),(0,0)"),
        7: (11, "(10,0),(15,0),(15,5),(15,-15),(9,3),(5,-5),(0,0)"),
        8: (17, "(0,10),(0,0),(5,5),(0,5),(5,10),(3,6),(3,9),(0,15)"),
        "9a": (25, "(0,1),(0,3),(1,0),(1,4),(2,2),(3,0),(3,4),(4,1),(4,3)"),
        "9b": (25, "(0,85),(0,-85),(-51,68),(51,-68),(-68,51),(68,-51),(-75,40),(75,-40),(0,0)"),
        10: (33, "(0,1105),(0,-1105),(-663,884),(663,-884),(-884,663),(884,-663),(-975,520),(975,-520),(-1020,425),(0,0)"),
        11: (41, "(0,1105),(0,-1105),(-663,884),(663,-884),(-884,663),(884,-663),(-975,520),(975,-520),(-1020,425),(1020,-425),(0,0)"),
        12: (51, "(0,40885),(0,-40885),(-24531,32708),(24531,-32708),(-32708,24531),(32708,-24531),(-36075,19240),(36075,-19240),(-37740,15725),(37740,-15725),(-38675,13260),(0,0)"),
        13: (61, "(0,40885),(0,-40885),(-24531,32708),(24531,-32708),(-32708,24531),(32708,-24531),(-36075,19240),(36075,-19240),(-37740,15725),(37740,-15725),(-38675,13260),(38675,-13260),(0,0)"),
        14: (73, "(0,204425),(0,-204425),(-122655,163540),(122655,-163540),(-163540,122655),(163540,-122655),(-180375,96200),(180375,-96200),(-188700,78625),(188700,-78625),(-193375,66300),(193375,-66300),(-196248,57239),(0,0)"),
        15: (85, "(0,204425),(0,-204425),(-122655,163540),(122655,-163540),(-163540,122655),(163540,-122655),(-180375,96200),(180375,-96200),(-188700,78625),(188700,-78625),(-193375,66300),(193375,-66300),(-196248,57239),(196248,-57239),(0,0)"),
        16: (99, "(0,204425),(0,-204425),(-122655,163540),(122655,-163540),(-163540,122655),(163540,-122655),(-180375,96200),(180375,-96200),(-188700,78625),(188700,-78625),(-193375,66300),(193375,-66300),(-196248,57239),(196248,-57239),(-198135,50320),(0,0)"),
    }

    def parse(s):
        s = s.replace(" ", "")
        return [tuple(int(t) for t in pair.split(",")) for pair in s[1:-1].split("),(")]

    ok_all = True
    print("== records from the task statement ==")
    for label, (claimed, s) in REC.items():
        pts = parse(s)
        r = analyse(pts)
        ok = (r["circles"] == claimed and not r["all_collinear"] and not r["all_concyclic"]
              and claimed <= f(r["n"]) + (0 if r["n"] >= 9 else 10))
        ok_all &= ok
        print(f"n={r['n']:2d} label={label!s:3} claimed={claimed:3d} mine={r['circles']:3d} f={f(r['n']):3d} "
              f"lines={r['lines']} circle_sizes={dict(Counter(r['circle_sizes']))} line_sizes={dict(Counter(r['line_sizes']))} "
              f"degenerate={r['all_collinear'] or r['all_concyclic']}  {'OK' if ok else 'MISMATCH'}")

    print("\n== search-found coordinates in search-local/runs/records.json ==")
    their = json.load(open(sys.argv[1])) if len(sys.argv) > 1 else {}
    for k, v in their.items():
        pts = [tuple(c) for c in v["coords"]]
        r = analyse(pts)
        ok = r["circles"] == v["circles"] and not r["all_collinear"] and not r["all_concyclic"]
        ok_all &= ok
        print(f"n={r['n']:2d} claimed={v['circles']:3d} mine={r['circles']:3d} f={f(r['n']):3d} lines={r['lines']} "
              f"circle_sizes={dict(Counter(r['circle_sizes']))} line_sizes={dict(Counter(r['line_sizes']))} "
              f"their_sizes={v['circle_sizes']}/{v['line_sizes']}  {'OK' if ok else 'MISMATCH'}")
    print("\nALL OK" if ok_all else "\nSOME MISMATCH")
