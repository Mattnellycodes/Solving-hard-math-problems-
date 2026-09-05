"""Fresh exact counter for Erdos #506 (audit v2, written independently).

Every triple of distinct points determines a unique generalised circle
    a (x^2+y^2) + b x + c y + d = 0
with (a,b,c,d) given by the 3x3 minors of the 3x4 matrix [x^2+y^2, x, y, 1]
(row per point).  a == 0  <=> the three points are collinear.  We normalise the
integer vector (a,b,c,d) by gcd and sign, so a generalised circle is identified
by a canonical integer key; the number of distinct keys with a != 0 is circles(P)
and with a == 0 is lines(P).  Integer-only arithmetic (inputs may be Fractions;
they are cleared to integers first).
"""
from fractions import Fraction
from itertools import combinations
from math import gcd, comb, lcm
from collections import defaultdict

def _to_int_points(points):
    pts = [(Fraction(x), Fraction(y)) for x, y in points]
    den = 1
    for x, y in pts:
        den = lcm(den, x.denominator, y.denominator)
    return [(int(x * den), int(y * den)) for x, y in pts]

def gcircle_key(p, q, r):
    rows = [(x * x + y * y, x, y, 1) for x, y in (p, q, r)]
    def minor(cols):
        (a1, a2, a3), (b1, b2, b3), (c1, c2, c3) = [tuple(row[c] for c in cols) for row in rows]
        return a1 * (b2 * c3 - b3 * c2) - a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1)
    a = minor((1, 2, 3)); b = -minor((0, 2, 3)); c = minor((0, 1, 3)); d = -minor((0, 1, 2))
    g = gcd(gcd(abs(a), abs(b)), gcd(abs(c), abs(d)))
    assert g != 0, "coincident points"
    a, b, c, d = a // g, b // g, c // g, d // g
    for v in (a, b, c, d):
        if v != 0:
            if v < 0:
                a, b, c, d = -a, -b, -c, -d
            break
    return (a, b, c, d)

def on_gcircle(key, p):
    a, b, c, d = key
    x, y = p
    return a * (x * x + y * y) + b * x + c * y + d == 0

def analyse(points):
    P = _to_int_points(points)
    n = len(P)
    assert len(set(P)) == n, "repeated point"
    blocks = defaultdict(set)
    for i, j, k in combinations(range(n), 3):
        blocks[gcircle_key(P[i], P[j], P[k])].update((i, j, k))
    # sanity: membership recomputed independently for every block
    for key, members in blocks.items():
        mem2 = {i for i in range(n) if on_gcircle(key, P[i])}
        assert mem2 == members, "inconsistent block membership"
    assert sum(comb(len(m), 3) for m in blocks.values()) == comb(n, 3)
    circ = {k: m for k, m in blocks.items() if k[0] != 0}
    lin = {k: m for k, m in blocks.items() if k[0] == 0}
    csz = sorted((len(m) for m in circ.values()), reverse=True)
    lsz = sorted((len(m) for m in lin.values()), reverse=True)
    degenerate = (bool(csz) and csz[0] == n) or (bool(lsz) and lsz[0] == n)
    D = sum(comb(k, 3) - 1 for k in csz + lsz if k >= 4)
    return dict(n=n, circles=len(circ), lines=len(lin), circle_sizes=csz, line_sizes=lsz,
                degenerate=degenerate, D=D, check_formula=comb(n, 3) - D - len(lin))

def f(n):
    return comb(n - 1, 2) + 1 - (n - 1) // 2
