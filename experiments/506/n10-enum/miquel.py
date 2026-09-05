#!/usr/bin/env python3
"""miquel.py -- Miquel closure test for abstract Moebius block structures, plus a numerical sanity check of
Miquel's theorem itself.

Miquel's theorem (real inversive plane): let p_v (v in {0,1}^3) be eight distinct points.  If five of the six
"faces" {v : v_i = c} are concyclic (each lies on a circle or line) then so is the sixth.  Proof sketch with
complex cross-ratios: a quadruple is concyclic iff its cross-ratio is real, and the product of suitably
arranged face cross-ratios is an identity, so the sixth is real if the other five are.  We verify the
theorem numerically below (random configurations with five faces forced concyclic, including the
degenerate case where two adjacent faces lie on the same circle).

Consequence for an abstract structure (every triple in exactly one block): if five faces of a cube on eight
points of P are each contained in a block, the sixth face must be contained in a block too.  A structure
violating this is not realisable.
Usage: python3 miquel.py results.json
"""
import sys, json, itertools, random, cmath


def cube_faces():
    verts = list(itertools.product((0, 1), repeat=3))
    faces = []
    for i in range(3):
        for c in (0, 1):
            faces.append([verts.index(v) for v in verts if v[i] == c])
    return verts, faces


VERTS, FACES = cube_faces()


def cube_labellings(points):
    """all injective labellings of the cube vertices by 8 points, modulo the cube's symmetry group (48)."""
    seen = set()
    out = []
    for perm in itertools.permutations(points):
        key = frozenset(frozenset(perm[i] for i in f) for f in FACES)
        if key in seen:
            continue
        seen.add(key)
        out.append([frozenset(perm[i] for i in f) for f in FACES])
    return out


def miquel_violations(n, blocks):
    blocks = [frozenset(b) for b in blocks]
    viol = []
    cache = {}
    for pts in itertools.combinations(range(n), 8):
        for faces in cube_labellings(pts):
            inblock = [any(f <= b for b in blocks) for f in faces]
            if sum(inblock) == 5:
                viol.append([sorted(f) for f in faces])
    return viol


def numeric_check(trials=200):
    """random Miquel configurations: five faces concyclic => sixth concyclic (|Im cross-ratio| ~ 0)."""
    def circle_through(a, b, c):
        # returns function t -> point on circle through a, b, c (parametrised by Moebius image of real t)
        # point on circle: w with cross ratio (w, a; b, c) real: w = (a*t*(b - c) - b*(a - c)) / (t*(b - c) - (a - c))
        return lambda t: (a * t * (b - c) - b * (a - c)) / (t * (b - c) - (a - c))

    def cr(a, b, c, d):
        return ((a - c) * (b - d)) / ((a - d) * (b - c))

    worst = 0.0
    worst_deg = 0.0
    rnd = lambda: complex(random.uniform(-2, 2), random.uniform(-2, 2))
    for k in range(trials):
        degenerate = (k % 2 == 1)
        # vertices: v000 v001 v010 v011 v100 v101 v110 v111 indexed 0..7 as in VERTS
        p = [None] * 8
        p[0], p[1], p[2] = rnd(), rnd(), rnd()
        # face x=0: {0,1,2,3} concyclic
        p[3] = circle_through(p[0], p[1], p[2])(random.uniform(-3, 3))
        # face y=0: {0,1,4,5} concyclic
        if degenerate:
            p[4] = circle_through(p[0], p[1], p[2])(random.uniform(-3, 3))   # same circle as face x=0
        else:
            p[4] = rnd()
        p[5] = circle_through(p[0], p[1], p[4])(random.uniform(-3, 3))
        # face z=0: {0,2,4,6} concyclic
        p[6] = circle_through(p[0], p[2], p[4])(random.uniform(-3, 3))
        # face z=1: {1,3,5,7} and face y=1: {2,3,6,7}: p7 = second intersection of circles (1,3,5) and (2,3,6)
        # solve: cross ratio (p7,1;3,5) real and (p7,2;3,6) real; parametrise p7 on circle(1,3,5) by t and solve
        # for t with Im cr(p7, 2, 3, 6) = 0 -> a real rational equation; use numeric root finding on t.
        f = circle_through(p[1], p[3], p[5])
        g = lambda t: cr(f(t), p[2], p[3], p[6]).imag
        # scan for sign changes
        ts = [i / 50.0 for i in range(-500, 501)]
        roots = []
        for t0, t1 in zip(ts, ts[1:]):
            try:
                g0, g1 = g(t0), g(t1)
            except ZeroDivisionError:
                continue
            if g0 == 0 or g0 * g1 < 0:
                lo, hi = t0, t1
                for _ in range(60):
                    mid = (lo + hi) / 2
                    if g(lo) * g(mid) <= 0:
                        hi = mid
                    else:
                        lo = mid
                roots.append((lo + hi) / 2)
        cand = [f(t) for t in roots if abs(f(t) - p[3]) > 1e-6]
        if not cand:
            continue
        p[7] = cand[0]
        pts = p
        # check all six faces
        vals = []
        for face in FACES:
            a, b, c, d = [pts[i] for i in face]
            vals.append(abs(cr(a, b, c, d).imag))
        if max(vals[:5]) > 1e-6:
            continue
        if degenerate:
            worst_deg = max(worst_deg, vals[5])
        else:
            worst = max(worst, vals[5])
    return worst, worst_deg


if __name__ == "__main__":
    random.seed(1)
    w, wd = numeric_check()
    print(f"numerical check of Miquel's theorem: max |Im cr(6th face)| = {w:.2e} (generic), {wd:.2e} (two faces on one circle)")
    if len(sys.argv) > 1:
        data = json.load(open(sys.argv[1]))
        recs = data['results'] if isinstance(data, dict) else data
        for i, rec in enumerate(recs):
            v = miquel_violations(rec['n'], rec['blocks'])
            print(f"structure {i}: sizes={rec['sizes']} count_min={rec['count_min']}: Miquel violations = {len(v)}"
                  + (f"; e.g. faces {v[0]}" if v else "  (Miquel-closed)"))
