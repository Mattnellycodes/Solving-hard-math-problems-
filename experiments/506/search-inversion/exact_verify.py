"""Exact certification of circle counts.

count_circles_exact(points, K)  -- points with coordinates in a field K (python Fractions, or a sympy
                                   algebraic number field QQ.algebraic_field(...)); returns the number
                                   of circles through >= 3 points, the number of collinear triples, the
                                   block sizes, and a degeneracy flag.  All arithmetic is exact.
invert_exact(points, O, K)      -- exact inversion (unit radius) about O.

Independent of ../circles_exact.py (different code path, generic field), so that records can be
double-checked with both.
"""
from fractions import Fraction as Fr
from itertools import combinations
from collections import defaultdict
import sympy
from sympy import QQ, sqrt, Rational, nsimplify


class RationalField:
    """Adapter so that plain Fractions look like a sympy domain."""
    def from_sympy(self, e):
        e = sympy.nsimplify(e)
        if not e.is_Rational:
            raise ValueError(f"not rational: {e}")
        return Fr(int(e.p), int(e.q))
    def to_sympy(self, a):
        return Rational(a.numerator, a.denominator)
    zero = Fr(0)
    def __call__(self, a):
        return Fr(a)


def make_field(radicands=()):
    """K = Q(sqrt(r1), sqrt(r2), ...) as a sympy algebraic field (or the rationals)."""
    radicands = [r for r in radicands if r not in (0, 1)]
    if not radicands:
        return RationalField()
    gens = [sqrt(Rational(r)) for r in radicands]
    return QQ.algebraic_field(*gens)


def conv(K, e):
    """Convert a sympy expression / Fraction / int to an element of K."""
    if isinstance(K, RationalField):
        if isinstance(e, Fr):
            return e
        return K.from_sympy(sympy.sympify(e))
    if isinstance(e, Fr):
        e = Rational(e.numerator, e.denominator)
    return K.from_sympy(sympy.sympify(e))


def zero_of(K):
    return Fr(0) if isinstance(K, RationalField) else K.zero


def circle_key(a, b, c, K):
    """Exact circle through a,b,c as a hashable key (cx, cy, r2), or ('L', A, B, C) for a line,
    with the line normalised so that the key is canonical."""
    ax, ay = a; bx, by = b; cx, cy = c
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if d == zero_of(K):
        # line through a,b: A x + B y = C with (A,B,C) normalised (first nonzero of A,B is 1)
        A = by - ay; B = ax - bx; C = A * ax + B * ay
        if A != zero_of(K):
            B = B / A; C = C / A; A = A / A
        else:
            C = C / B; B = B / B
        return ('L', A, B, C)
    a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx * cx + cy * cy
    ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
    uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
    r2 = (ax - ux) ** 2 + (ay - uy) ** 2
    return ('C', ux, uy, r2)


def count_circles_exact(points, K):
    """points: list of (x,y) with x,y elements of K.  Returns dict."""
    n = len(points)
    assert len(set(points)) == n, "duplicate points"
    blocks = defaultdict(set)
    for i, j, k in combinations(range(n), 3):
        key = circle_key(points[i], points[j], points[k], K)
        blocks[key].update((i, j, k))
    circles = {k: v for k, v in blocks.items() if k[0] == 'C'}
    lines = {k: v for k, v in blocks.items() if k[0] == 'L'}
    degenerate = any(len(v) == n for v in blocks.values())
    sizes = sorted(len(v) for v in circles.values())
    return dict(n=n, circles=len(circles), lines=len(lines), blocks=len(blocks),
                circle_sizes=sizes, line_sizes=sorted(len(v) for v in lines.values()),
                degenerate=degenerate)


def invert_exact(points, O, K):
    """Unit-radius inversion about O (an element pair of K); O must not be one of the points."""
    ox, oy = O
    out = []
    for (x, y) in points:
        dx = x - ox; dy = y - oy
        r2 = dx * dx + dy * dy
        if r2 == zero_of(K):
            raise ValueError("inversion centre coincides with a point")
        out.append((ox + dx / r2, oy + dy / r2))
    return out


def certify(points_sympy, O_sympy=None, radicands=()):
    """points_sympy: list of (x,y) sympy expressions (or 'inf' allowed only if O is given);
    O_sympy: optional inversion centre.  Returns the exact count dict of the (inverted) set and the
    exact coordinates (as sympy expressions) of the certified planar set."""
    K = make_field(radicands)
    pts = []
    has_inf = False
    for p in points_sympy:
        if p == 'inf':
            has_inf = True; continue
        pts.append((conv(K, p[0]), conv(K, p[1])))
    if O_sympy is not None:
        O = (conv(K, O_sympy[0]), conv(K, O_sympy[1]))
        inv = invert_exact(pts, O, K)
        if has_inf:
            inv.append(O)      # infinity maps to the centre
        pts = inv
    elif has_inf:
        raise ValueError("'inf' in the set but no inversion centre")
    res = count_circles_exact(pts, K)
    if isinstance(K, RationalField):
        coords = [(Rational(x.numerator, x.denominator), Rational(y.numerator, y.denominator)) for x, y in pts]
    else:
        coords = [(K.to_sympy(x), K.to_sympy(y)) for x, y in pts]
    res['coords'] = coords
    return res


if __name__ == "__main__":
    # 1. rational record n=8
    pts = [(0, 0), (0, 5), (0, 10), (0, 15), (3, 6), (3, 9), (5, 5), (5, 10)]
    print("n=8 record:", {k: v for k, v in certify(pts).items() if k != 'coords'})
    # 2. orthocentric 7 + inf inverted in O=(0,1) -> must give 17
    pts7 = [(0, 3), (1, 0), (3, 0), (0, -1), (Rational(6, 5), Rational(-3, 5)), (2, 1), (0, 0), 'inf']
    r = certify(pts7, O_sympy=(0, 1))
    print("orthocentric+inf inverted at (0,1):", {k: v for k, v in r.items() if k != 'coords'})
    print("   coords:", r['coords'])
    # 3. an algebraic example: two concentric squares, ratio 4+sqrt(15) (symmetric agent: 17)
    rho = 4 + sqrt(15)
    sq = [(1, 0), (0, 1), (-1, 0), (0, -1), (rho, 0), (0, rho), (-rho, 0), (0, -rho)]
    r = certify(sq, radicands=(15,))
    print("two squares ratio 4+sqrt15 (no inversion):", {k: v for k, v in r.items() if k != 'coords'})
    # 4. antipodal configuration n=10 (rational points on the unit circle)
    from fractions import Fraction
    def formula_config(n):
        pts = [(0, 0)]
        t = 1
        while len(pts) < n:
            x = Rational(1 - t * t, 1 + t * t); y = Rational(2 * t, 1 + t * t)
            for p in [(x, y), (-x, -y)]:
                if len(pts) < n and p not in pts:
                    pts.append(p)
            t += 1
        return pts
    for n in range(9, 17):
        r = certify(formula_config(n))
        print(f"antipodal n={n}: circles={r['circles']} formula={(n-1)*(n-2)//2+1-(n-1)//2}")
