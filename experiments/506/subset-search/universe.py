"""Universes of candidate points for the subset search (Erdos #506).

A universe is a finite point set U (given exactly and as floats) together with the point at
infinity (index N = len(U)).  Its *blocks* are the circles and lines through >= 3 points of U;
lines contain the point at infinity.  Because every triple of U lies in exactly one block, the
Moebius block structure of any subset S of U u {inf} is {B : |B & S| >= 3}, so the circle count of the
best inversion of S is  |B(S)| - max_{O in (U u {inf}) \ S} deg_S(O)  (a candidate-restricted upper
bound on the true optimum, which is re-evaluated later with all block intersections).
"""
import math, itertools, json, os
from fractions import Fraction as Fr
from math import gcd
import numpy as np
import sympy as sp

# ---------------------------------------------------------------- exact block keys (rational)

def _norm_rat_pair(x, y):
    return (x, y)

def blocks_rational(pts):
    """pts: list of (Fraction, Fraction).  Returns list of blocks (sorted tuples of indices),
    lines get the extra index N (= infinity)."""
    N = len(pts)
    # scale to integers
    den = 1
    for x, y in pts:
        den = den * x.denominator // gcd(den, x.denominator)
        den = den * y.denominator // gcd(den, y.denominator)
    P = [(int(x * den), int(y * den)) for x, y in pts]
    assert len(set(P)) == N, "duplicate points"
    blocks = {}
    for i in range(N):
        xi, yi = P[i]
        for j in range(i + 1, N):
            xj, yj = P[j]
            dx1, dy1 = xj - xi, yj - yi
            s1 = xj * xj - xi * xi + yj * yj - yi * yi
            for k in range(j + 1, N):
                xk, yk = P[k]
                dx2, dy2 = xk - xi, yk - yi
                det = dx1 * dy2 - dx2 * dy1
                if det == 0:
                    A, B = dy1, -dx1
                    C = A * xi + B * yi
                    g = gcd(gcd(abs(A), abs(B)), abs(C))
                    A //= g; B //= g; C //= g
                    if A < 0 or (A == 0 and B < 0): A, B, C = -A, -B, -C
                    key = ('L', A, B, C)
                else:
                    s2 = xk * xk - xi * xi + yk * yk - yi * yi
                    # centre = (s1*dy2 - s2*dy1, dx1*s2 - dx2*s1) / (2 det)
                    un = s1 * dy2 - s2 * dy1; vn = dx1 * s2 - dx2 * s1; d = 2 * det
                    g = gcd(gcd(abs(un), abs(vn)), abs(d))
                    un //= g; vn //= g; d //= g
                    if d < 0: un, vn, d = -un, -vn, -d
                    r2n = (xi * d - un) ** 2 + (yi * d - vn) ** 2
                    key = ('C', un, vn, d, r2n)
                s = blocks.get(key)
                if s is None:
                    blocks[key] = {i, j, k}
                else:
                    s.update((i, j, k))
    out = []
    for key, s in blocks.items():
        if key[0] == 'L': s.add(N)
        out.append(tuple(sorted(s)))
    return add_two_point_lines(out, N)

def add_two_point_lines(out, N):
    """lines through exactly two finite points are 3-point blocks of U u {inf}: add them."""
    covered = set()
    for b in out:
        if b[-1] == N:
            fin = b[:-1]
            for i in range(len(fin)):
                for j in range(i + 1, len(fin)):
                    covered.add((fin[i], fin[j]))
    for i in range(N):
        for j in range(i + 1, N):
            if (i, j) not in covered:
                out.append((i, j, N))
    return out

# ---------------------------------------------------------------- float block keys (algebraic)

def blocks_float(P, tol=1e-7):
    """P: (N,2) float array. Groups triples by circumcircle/line with tolerance clustering."""
    N = len(P)
    keys = {}
    trip_key = []
    # first pass: rounded keys with clustering via union-find over near keys
    raw = {}
    order = []
    for i, j, k in itertools.combinations(range(N), 3):
        ax, ay = P[i]; bx, by = P[j]; cx, cy = P[k]
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        scale = max(abs(bx-ax), abs(by-ay), abs(cx-ax), abs(cy-ay), 1e-300)
        if abs(d) < 1e-9 * scale * scale * 4 + 1e-12:
            A = by - ay; B = ax - bx; C = A * ax + B * ay
            nrm = math.hypot(A, B); A /= nrm; B /= nrm; C /= nrm
            if A < -1e-9 or (abs(A) < 1e-9 and B < 0): A, B, C = -A, -B, -C
            key = (0.0, A, B, C)
        else:
            a2 = ax*ax+ay*ay; b2 = bx*bx+by*by; c2 = cx*cx+cy*cy
            ux = (a2*(by-cy)+b2*(cy-ay)+c2*(ay-by))/d
            uy = (a2*(cx-bx)+b2*(ax-cx)+c2*(bx-ax))/d
            r = math.hypot(ax-ux, ay-uy)
            key = (1.0, ux, uy, r)
        raw.setdefault(tuple(round(v / tol) for v in key), []).append((i, j, k))
    # merge neighbouring cells
    parent = {}
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    cells = list(raw.keys())
    for c in cells: parent[c] = c
    for c in cells:
        for d0 in (0,):
            for d1 in (-1, 0, 1):
                for d2 in (-1, 0, 1):
                    for d3 in (-1, 0, 1):
                        nb = (c[0], c[1]+d1, c[2]+d2, c[3]+d3)
                        if nb != c and nb in parent:
                            ra, rb = find(c), find(nb)
                            if ra != rb: parent[ra] = rb
    groups = {}
    for c in cells:
        groups.setdefault(find(c), set()).update(itertools.chain.from_iterable(raw[c]))
    out = []
    for c, s in groups.items():
        if c[0] == 0: s.add(N)
        out.append(tuple(sorted(s)))
    return add_two_point_lines(out, N)

# ---------------------------------------------------------------- universe class

class Universe:
    def __init__(self, name, exact, floats=None, blocks=None):
        """exact: list of (x,y) with Fraction or sympy entries (strings allowed)."""
        self.name = name
        self.exact = exact
        if floats is None:
            floats = [(float(x), float(y)) for x, y in exact]
        self.P = np.array(floats, dtype=float)
        self.N = len(self.P)
        self.rational = all(isinstance(x, (int, Fr)) and isinstance(y, (int, Fr)) for x, y in exact)
        if blocks is None:
            if self.rational:
                blocks = blocks_rational([(Fr(x), Fr(y)) for x, y in exact])
            else:
                blocks = blocks_float(self.P)
        self.blocks = blocks
        self.build_csr()

    def build_csr(self):
        N1 = self.N + 1
        M = len(self.blocks)
        bsize = np.array([len(b) for b in self.blocks], dtype=np.int32)
        bptr = np.zeros(M + 1, dtype=np.int64); bptr[1:] = np.cumsum(bsize)
        bmem = np.concatenate([np.array(b, dtype=np.int32) for b in self.blocks]) if M else np.zeros(0, np.int32)
        pl = [[] for _ in range(N1)]
        for bi, b in enumerate(self.blocks):
            for p in b: pl[p].append(bi)
        pptr = np.zeros(N1 + 1, dtype=np.int64); pptr[1:] = np.cumsum([len(x) for x in pl])
        pblk = np.concatenate([np.array(x, dtype=np.int32) for x in pl]) if N1 else np.zeros(0, np.int32)
        self.bptr, self.bmem, self.pptr, self.pblk = bptr, bmem, pptr, pblk

    def exact_str(self, i):
        if i == self.N: return "inf"
        x, y = self.exact[i]
        return f"({x}, {y})"

    def summary(self):
        from collections import Counter
        sizes = Counter(len(b) for b in self.blocks)
        lines = sum(1 for b in self.blocks if b[-1] == self.N)
        return f"{self.name}: N={self.N} points, {len(self.blocks)} blocks ({lines} lines), sizes {dict(sorted(sizes.items()))}"

# ---------------------------------------------------------------- constructors

def grid(k, l=None, step=1):
    l = l or k
    pts = [(Fr(i, step), Fr(j, step)) for i in range(k) for j in range(l)]
    return Universe(f"grid{k}x{l}" + (f"_step1/{step}" if step != 1 else ""), pts)

def grid_half(k):
    """integer points of [0,k-1]^2 plus the centres of the unit squares (half-integer points)."""
    pts = [(Fr(i), Fr(j)) for i in range(k) for j in range(k)]
    pts += [(Fr(2*i+1, 2), Fr(2*j+1, 2)) for i in range(k-1) for j in range(k-1)]
    return Universe(f"grid{k}x{k}+halfcentres", pts)

def lattice_circles(Ns, centre=True):
    """lattice points on the concentric circles x^2+y^2 = N for N in Ns (+ origin)."""
    pts = []
    for Nn in Ns:
        r = int(math.isqrt(Nn))
        for x in range(-r, r + 1):
            y2 = Nn - x * x
            y = int(math.isqrt(y2))
            if y * y == y2:
                pts.append((Fr(x), Fr(y)))
                if y: pts.append((Fr(x), Fr(-y)))
    if centre: pts.append((Fr(0), Fr(0)))
    pts = sorted(set(pts))
    return Universe(f"lattice_circles{list(Ns)}" + ("+c" if centre else ""), pts)

def polygons(m, radii, rots=(0, 1), centre=True, name=None):
    """concentric regular m-gons with the given radii (sympy exprs), rotations 0 and pi/m."""
    ex = []; fl = []
    for r in radii:
        r = sp.nsimplify(r)
        for rot in rots:
            for k in range(m):
                ang = 2 * sp.pi * k / m + rot * sp.pi / m
                x = r * sp.cos(ang); y = r * sp.sin(ang)
                ex.append((x, y)); fl.append((float(x), float(y)))
    if centre:
        ex.append((sp.Integer(0), sp.Integer(0))); fl.append((0.0, 0.0))
    # dedupe
    keep = []; out_ex = []; out_fl = []
    for e, f in zip(ex, fl):
        if all(abs(f[0]-g[0]) > 1e-9 or abs(f[1]-g[1]) > 1e-9 for g in out_fl):
            out_ex.append(e); out_fl.append(f)
    return Universe(name or f"poly{m}_r{[str(r) for r in radii]}_rot{list(rots)}" + ("+c" if centre else ""), out_ex, out_fl)

def circle_plus_grid(m, R=1, grid_range=2, extra_centre=True):
    """regular m-gon of radius R (sympy) plus the integer points of [-g,g]^2 (as sympy)."""
    ex = []; fl = []
    for k in range(m):
        ang = 2 * sp.pi * k / m
        x = sp.nsimplify(R) * sp.cos(ang); y = sp.nsimplify(R) * sp.sin(ang)
        ex.append((x, y)); fl.append((float(x), float(y)))
    for i in range(-grid_range, grid_range + 1):
        for j in range(-grid_range, grid_range + 1):
            f = (float(i), float(j))
            if all(abs(f[0]-g[0]) > 1e-9 or abs(f[1]-g[1]) > 1e-9 for g in fl):
                ex.append((sp.Integer(i), sp.Integer(j))); fl.append(f)
    return Universe(f"{m}gon_R{R}+grid[-{grid_range},{grid_range}]^2", ex, fl)

def augment(U, min_deg=4, max_add=60, name=None):
    """Add to U the pairwise intersection points of its blocks that lie on >= min_deg blocks
    (float geometry; exact coordinates are recovered later for records)."""
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
    from explore_families import circ_line_intersections, EPS, key
    P = U.P; N = U.N
    geo = []
    for b in U.blocks:
        idx = [p for p in b if p != N][:3]
        a, bb, c = P[idx[0]], P[idx[1]], P[idx[2]]
        ax, ay = a; bx, by = bb; cx, cy = c
        d = 2 * (ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
        if b[-1] == N:
            A = by-ay; B = ax-bx; C = A*ax+B*ay
            nrm = math.hypot(A, B); A/=nrm; B/=nrm; C/=nrm
            if A < -1e-12 or (abs(A) < 1e-12 and B < 0): A, B, C = -A, -B, -C
            geo.append(('L', key(A), key(B), key(C)))
        else:
            a2=ax*ax+ay*ay; b2=bx*bx+by*by; c2=cx*cx+cy*cy
            ux=(a2*(by-cy)+b2*(cy-ay)+c2*(ay-by))/d; uy=(a2*(cx-bx)+b2*(ax-cx)+c2*(bx-ax))/d
            geo.append(('C', key(ux), key(uy), key(math.hypot(ax-ux, ay-uy))))
    ptkeys = set((round(x/1e-5), round(y/1e-5)) for x, y in P)
    seen = {}
    for i in range(len(geo)):
        for j in range(i+1, len(geo)):
            for (x, y) in circ_line_intersections(geo[i], geo[j]):
                kk = (round(x/1e-5), round(y/1e-5))
                if kk in ptkeys: continue
                s = seen.setdefault(kk, set()); s.add(i); s.add(j)
    cand = sorted(((len(s), kk) for kk, s in seen.items() if len(s) >= min_deg), reverse=True)[:max_add]
    newfl = list(map(tuple, P)) + [(kk[0]*1e-5, kk[1]*1e-5) for _, kk in cand]
    # refine the float coordinates of new points by intersecting two of their blocks exactly-ish:
    # (we keep 1e-5-rounded values; records are re-derived exactly later)
    newex = list(U.exact) + [(sp.Float(x, 15), sp.Float(y, 15)) for x, y in newfl[N:]]
    return Universe(name or (U.name + f"+aug{len(cand)}"), newex, newfl)

if __name__ == "__main__":
    import time
    t = time.time()
    U = grid(5); print(U.summary(), time.time() - t)
    U = polygons(6, [1, 2]); print(U.summary(), time.time() - t)
