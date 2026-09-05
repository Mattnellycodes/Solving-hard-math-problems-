"""Adjudicator's own exact recount of the record configurations (written from scratch).

Representation deliberately different from every other counter in this experiment: no circle
equations, no centres/radii.  A block is identified purely by its POINT SET, computed from exact
integer 3x3 (collinearity) and 4x4 (concyclicity) determinants:
    block(i,j,k) = {q : q is concyclic with p_i,p_j,p_k}      if the triple is not collinear
    block(i,j,k) = {q : q is collinear with p_i,p_j}          if the triple is collinear
Distinct point sets = distinct blocks.  circles = #non-line blocks.  We then verify the Moebius
identity circles = C(n,3) - D - ell, the intersection axioms, and non-degeneracy.
"""
from itertools import combinations
from math import comb
from fractions import Fraction as Fr


def det(M):
    """Exact determinant by Laplace expansion (small matrices only)."""
    k = len(M)
    if k == 1:
        return M[0][0]
    s = 0
    for j in range(k):
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        s += (-1) ** j * M[0][j] * det(minor)
    return s


def collinear(p, q, r):
    return det([[p[0], p[1], 1], [q[0], q[1], 1], [r[0], r[1], 1]]) == 0


def concyclic(p, q, r, s):
    return det([[t[0] * t[0] + t[1] * t[1], t[0], t[1], 1] for t in (p, q, r, s)]) == 0


def analyse(points):
    P = [(Fr(x), Fr(y)) for x, y in points]
    n = len(P)
    assert len(set(P)) == n, "repeated point"
    blocks = {}  # frozenset -> 'L' or 'C'
    for i, j, k in combinations(range(n), 3):
        if collinear(P[i], P[j], P[k]):
            S = frozenset(q for q in range(n) if collinear(P[i], P[j], P[q]))
            kind = 'L'
        else:
            S = frozenset({i, j, k} | {q for q in range(n) if q not in (i, j, k)
                                        and concyclic(P[i], P[j], P[k], P[q])})
            kind = 'C'
            # sanity: no three points of a circle block are collinear
            for a, b, c in combinations(sorted(S), 3):
                assert not collinear(P[a], P[b], P[c]), "circle block contains a collinear triple"
        if S in blocks:
            assert blocks[S] == kind
        blocks[S] = kind
    # every triple lies in exactly one block
    cover = {t: 0 for t in combinations(range(n), 3)}
    for S in blocks:
        for t in combinations(sorted(S), 3):
            cover[t] += 1
    assert all(v == 1 for v in cover.values()), "triple covered != 1 times"
    # intersection axioms
    Bl = list(blocks)
    for A, B in combinations(Bl, 2):
        assert len(A & B) <= 2
        if blocks[A] == 'L' and blocks[B] == 'L':
            assert len(A & B) <= 1
    circles = [S for S in Bl if blocks[S] == 'C']
    lines = [S for S in Bl if blocks[S] == 'L']
    D = sum(comb(len(S), 3) - 1 for S in Bl if len(S) >= 4)
    ell = len(lines)
    degenerate = any(len(S) == n for S in Bl)
    ident = (len(circles) == comb(n, 3) - D - ell)
    return dict(n=n, circles=len(circles), lines=ell, D=D, identity=ident, degenerate=degenerate,
                circle_sizes=sorted(len(S) for S in circles), line_sizes=sorted(len(S) for S in lines),
                big_blocks=sorted(sorted(S) for S in Bl if len(S) >= 4),
                line_blocks=sorted(sorted(S) for S in lines))


def f(n):
    return comb(n - 1, 2) + 1 - (n - 1) // 2


CONF = {
    'n=6 report': [(0, 0), (20, 0), (5, 15), (5, 0), (2, 6), (10, 10)],
    'n=7 report': [(0, 0), (20, 0), (5, 15), (5, 0), (2, 6), (10, 10), (5, 5)],
    'n=8 report': [(0, 0), (0, 5), (0, 10), (0, 15), (3, 6), (3, 9), (5, 5), (5, 10)],
    'n=8 two squares (control, 18)': [(1, 0), (0, 1), (-1, 0), (0, -1), (2, 0), (0, 2), (-2, 0), (0, -2)],
    'n=9 antipodal': [(0, 0), (5, 0), (-5, 0), (0, 5), (0, -5), (3, 4), (-3, -4), (3, -4), (-3, 4)],
    'n=9 8 concyclic + point on 4 non-diameter chords': None,
}

if __name__ == '__main__':
    # build the second n=9 witness: 8 rational points on the unit circle in pairs collinear with (1/3,0)
    # chord through (1/3,0) and unit-circle point t -> second intersection computed exactly
    def second(pt, q):
        # line through q=(qx,qy) and pt on unit circle; other intersection with x^2+y^2=1
        (x1, y1), (qx, qy) = pt, q
        dx, dy = x1 - qx, y1 - qy
        # param q + s*(d): |q+s d|^2 = 1 -> s^2 |d|^2 + 2 s (q.d) + |q|^2 - 1 = 0 ; s=1 is a root
        a = dx * dx + dy * dy
        b = 2 * (qx * dx + qy * dy)
        s2 = -(b) / a - 1  # sum of roots = -b/a, one root is 1
        return (qx + s2 * dx, qy + s2 * dy)
    q = (Fr(1, 3), Fr(0))
    base = [(Fr(3, 5), Fr(4, 5)), (Fr(-3, 5), Fr(4, 5)), (Fr(5, 13), Fr(12, 13)), (Fr(0), Fr(1))]
    pts9 = [q]
    for b in base:
        pts9.append(b)
        pts9.append(second(b, q))
    CONF['n=9 8 concyclic + point on 4 non-diameter chords'] = pts9
    for name, pts in CONF.items():
        r = analyse(pts)
        print(f"{name}: n={r['n']} circles={r['circles']} f(n)={f(r['n'])} lines={r['lines']} D={r['D']} "
              f"identity={r['identity']} degenerate={r['degenerate']}")
        print(f"    circle sizes {r['circle_sizes']}  line sizes {r['line_sizes']}")
        print(f"    blocks of size>=4: {r['big_blocks']}")
        print(f"    lines: {r['line_blocks']}")
