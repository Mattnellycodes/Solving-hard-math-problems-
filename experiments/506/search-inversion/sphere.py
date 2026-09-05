"""Sphere (Moebius) model for Erdos #506 construction search.

A Moebius configuration is a finite set of points on the unit sphere S^2.  A *block* is a plane
containing >= 3 of the points (its intersection with the sphere is the circle through them).
Stereographic projection from a sphere point O (not in the set) maps the blocks through O to
straight lines and all other blocks to circles, so

    circles(projection from O) = |B| - deg(O),   deg(O) = number of blocks through O.

The best planar picture of a Moebius set is therefore obtained from a sphere point O lying on as
many blocks as possible; such points are intersections of two block circles, so all candidates are
found by intersecting pairs of block planes with the sphere.

Everything here is floating point (numpy).  Records are re-verified exactly elsewhere
(exact_verify.py).
"""
import numpy as np, itertools, math
from collections import defaultdict

TOL_ANGLE = 1e-7      # clustering tolerance for plane angles (radians)
TOL_PLANE = 1e-7      # |n.x - d| tolerance for point-on-plane tests
TOL_POINT = 1e-6      # merging tolerance for candidate points


def formula(n):
    return (n - 1) * (n - 2) // 2 + 1 - (n - 1) // 2


# ---------------------------------------------------------------- coordinates
def plane_to_sphere(pts):
    """Inverse stereographic projection from the north pole (0,0,1).  'inf' -> north pole."""
    out = []
    for p in pts:
        if p == 'inf' or p is None:
            out.append((0.0, 0.0, 1.0))
        else:
            x, y = float(p[0]), float(p[1])
            r2 = x * x + y * y
            out.append((2 * x / (r2 + 1), 2 * y / (r2 + 1), (r2 - 1) / (r2 + 1)))
    return np.array(out, dtype=float)


def rotation_taking_to_north(o):
    """Orthogonal matrix R with R o = (0,0,1)."""
    o = np.asarray(o, float); o = o / np.linalg.norm(o)
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(o, z); s = np.linalg.norm(v); c = float(o @ z)
    if s < 1e-14:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))


def sphere_to_plane(V, o):
    """Stereographic projection from the sphere point o (o must not be one of the points)."""
    R = rotation_taking_to_north(o)
    W = np.asarray(V, float) @ R.T
    out = []
    for w in W:
        den = 1.0 - w[2]
        if abs(den) < 1e-12:
            out.append(None)
        else:
            out.append((w[0] / den, w[1] / den))
    return out


def normalise_rows(V):
    V = np.asarray(V, float)
    return V / np.linalg.norm(V, axis=1)[:, None]


# ---------------------------------------------------------------- blocks
def _cluster_1d(vals, tol, period=None):
    """Cluster sorted values; returns list of index lists.  Handles a periodic coordinate."""
    idx = np.argsort(vals)
    groups = []
    cur = [idx[0]]
    for a, b in zip(idx[:-1], idx[1:]):
        if vals[b] - vals[a] <= tol:
            cur.append(b)
        else:
            groups.append(cur); cur = [b]
    groups.append(cur)
    if period is not None and len(groups) > 1:
        if (vals[groups[0][0]] + period) - vals[groups[-1][-1]] <= tol:
            groups[0] = groups[-1] + groups[0]; groups.pop()
    return groups


def find_blocks(V, tol=TOL_ANGLE):
    """Return list of blocks (frozensets of point indices) of the sphere point set V (N x 3).
    Robust: for every pair (i,j) the planes through the chord ij are parametrised by an angle mod pi
    and the third points are clustered by that angle."""
    V = np.asarray(V, float); N = len(V)
    blocks = set()
    for i in range(N):
        for j in range(i + 1, N):
            d = V[j] - V[i]; dn = np.linalg.norm(d)
            if dn < 1e-12:
                raise ValueError(f"duplicate points {i},{j}")
            d = d / dn
            # orthonormal basis of the plane perpendicular to d
            a = np.array([1.0, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1.0, 0])
            e1 = a - (a @ d) * d; e1 /= np.linalg.norm(e1); e2 = np.cross(d, e1)
            others = [k for k in range(N) if k != i and k != j]
            if not others:
                continue
            W = V[others] - V[i]
            nrm = np.cross(np.broadcast_to(d, W.shape), W)   # normal of plane (i,j,k)
            ln = np.linalg.norm(nrm, axis=1)
            good = ln > 1e-10  # k not on the chord line (impossible on a sphere unless duplicate)
            ang = np.arctan2(nrm[:, :] @ e2, nrm[:, :] @ e1)  # in (-pi, pi]
            ang = np.mod(ang, np.pi)                         # plane is unoriented: mod pi
            vals = ang
            groups = _cluster_1d(vals, tol, period=np.pi)
            for g in groups:
                if len(g) >= 1:
                    blk = frozenset([i, j] + [others[t] for t in g if good[t]])
                    if len(blk) >= 3:
                        blocks.add(blk)
    return sorted(blocks, key=lambda b: (-len(b), sorted(b)))


def plane_of(V, blk):
    """Unit normal n and offset d with n.x = d for the plane through the points of blk (least squares)."""
    P = V[sorted(blk)]
    c = P.mean(axis=0)
    U, S, Wt = np.linalg.svd(P - c)
    n = Wt[-1]
    if n[2] < 0 or (abs(n[2]) < 1e-12 and (n[1] < 0 or (abs(n[1]) < 1e-12 and n[0] < 0))):
        n = -n
    return n, float(n @ c)


def planes(V, blocks):
    Nn = np.zeros((len(blocks), 3)); D = np.zeros(len(blocks))
    for t, b in enumerate(blocks):
        Nn[t], D[t] = plane_of(V, b)
    return Nn, D


# ---------------------------------------------------------------- candidate centres
def candidate_points(Nn, D, max_pairs=3_000_000):
    """All intersection points of pairs of block planes with the unit sphere (M x 3 array)."""
    B = len(Nn)
    pts = []
    for i in range(B):
        n1, d1 = Nn[i], D[i]
        n2 = Nn[i + 1:]; d2 = D[i + 1:]
        if len(n2) == 0:
            continue
        # line of intersection: point x0 = a n1 + b n2, direction u = n1 x n2
        u = np.cross(np.broadcast_to(n1, n2.shape), n2)
        un = np.linalg.norm(u, axis=1)
        ok = un > 1e-9
        if not ok.any():
            continue
        n2 = n2[ok]; d2 = d2[ok]; u = u[ok] / un[ok][:, None]
        g = n2 @ n1                      # cos of angle between normals
        det = 1.0 - g * g
        a = (d1 - d2 * g) / det; b = (d2 - d1 * g) / det
        x0 = a[:, None] * n1[None, :] + b[:, None] * n2
        # |x0 + t u|^2 = 1 with x0 . u = 0  ->  t^2 = 1 - |x0|^2
        h2 = 1.0 - np.sum(x0 * x0, axis=1)
        okh = h2 > -1e-9
        h = np.sqrt(np.maximum(h2[okh], 0.0))
        x0 = x0[okh]; u = u[okh]
        pts.append(x0 + h[:, None] * u); pts.append(x0 - h[:, None] * u)
        if sum(len(p) for p in pts) > max_pairs:
            break
    if not pts:
        return np.zeros((0, 3))
    return np.vstack(pts)


def merge_points(P, tol=TOL_POINT, min_mult=1):
    """Merge near-duplicate points by grid hashing (two offset grids so that clusters split by a cell
    boundary are still found); returns (unique points, multiplicity).  Only cells with multiplicity
    >= min_mult are returned (a point on d >= 3 planes appears in C(d,2) >= 3 plane pairs)."""
    if len(P) == 0:
        return P, np.zeros(0, int)
    reps = []; mults = []
    M = int(4.0 / tol) + 3
    for shift in (0.0, 0.5):
        key = np.floor(P / tol + shift).astype(np.int64)
        h = (key[:, 0] + M // 2) * M * M + (key[:, 1] + M // 2) * M + (key[:, 2] + M // 2)
        order = np.argsort(h, kind='stable')
        hs = h[order]
        starts = np.r_[0, np.flatnonzero(np.diff(hs)) + 1]
        ends = np.r_[starts[1:], len(hs)]
        mult = ends - starts
        sel = mult >= min_mult
        reps.append(P[order[starts[sel]]]); mults.append(mult[sel])
    R = np.vstack(reps); Mu = np.concatenate(mults)
    if len(R) == 0:
        return R, Mu
    # merge representatives closer than tol (sorted sweep on x)
    keep = np.ones(len(R), bool)
    order = np.argsort(R[:, 0])
    Rs = R[order]
    for a in range(len(Rs)):
        if not keep[order[a]]:
            continue
        b = a + 1
        while b < len(Rs) and Rs[b, 0] - Rs[a, 0] <= 2 * tol:
            if keep[order[b]] and np.all(np.abs(Rs[b] - Rs[a]) <= 2 * tol):
                keep[order[b]] = False
            b += 1
    R = R[keep]; Mu = Mu[keep]
    return normalise_rows(R), Mu


def degrees(Cands, Nn, D, tol=TOL_PLANE, chunk=20000):
    """Number of block planes through each candidate point (vectorised, chunked)."""
    out = np.zeros(len(Cands), int)
    for s in range(0, len(Cands), chunk):
        C = Cands[s:s + chunk]
        R = np.abs(C @ Nn.T - D[None, :]) < tol
        out[s:s + chunk] = R.sum(axis=1)
    return out


def blocks_through(x, Nn, D, tol=TOL_PLANE):
    return [t for t in range(len(Nn)) if abs(Nn[t] @ x - D[t]) < tol]


# ---------------------------------------------------------------- evaluation
class Evaluation:
    def __init__(self, V, blocks, Nn, D, cands, cdeg, pdeg):
        self.V = V; self.blocks = blocks; self.Nn = Nn; self.D = D
        self.cands = cands; self.cdeg = cdeg; self.pdeg = pdeg
        self.nb = len(blocks)
        self.best_deg = int(cdeg.max()) if len(cdeg) else 0
        self.best_count = self.nb - self.best_deg

    def best_centres(self, k=5):
        order = np.argsort(-self.cdeg)[:k]
        return [(int(self.cdeg[t]), self.cands[t]) for t in order]

    def summary(self):
        sizes = defaultdict(int)
        for b in self.blocks:
            sizes[len(b)] += 1
        return dict(n=len(self.V), nblocks=self.nb, sizes=dict(sorted(sizes.items())),
                    best_deg=self.best_deg, best_count=self.best_count,
                    formula=formula(len(self.V)),
                    derived=sorted(int(self.nb - d) for d in self.pdeg)[:3])


def evaluate(V, want_candidates=True, max_block_pairs=3_000_000):
    """Full evaluation of a sphere point set V (N x 3, unit vectors).
    Returns Evaluation with candidate centres (points of the sphere on >= 2 blocks, not in V)."""
    V = normalise_rows(V)
    N = len(V)
    blocks = find_blocks(V)
    if any(len(b) == N for b in blocks):
        return None  # degenerate: all points on one circle
    Nn, D = planes(V, blocks)
    pdeg = np.array([sum(1 for b in blocks if i in b) for i in range(N)])
    if not want_candidates or len(blocks) < 2:
        return Evaluation(V, blocks, Nn, D, np.zeros((0, 3)), np.zeros(0, int), pdeg)
    C = candidate_points(Nn, D, max_pairs=max_block_pairs)
    U, mult = merge_points(C, min_mult=1 if len(C) < 200000 else 2)
    # drop candidates coinciding with configuration points
    if len(U):
        dist = np.min(np.linalg.norm(U[:, None, :] - V[None, :, :], axis=2), axis=1)
        keep = dist > 1e-6
        U = U[keep]
    cdeg = degrees(U, Nn, D) if len(U) else np.zeros(0, int)
    return Evaluation(V, blocks, Nn, D, U, cdeg, pdeg)


def evaluate_planar(pts, **kw):
    """pts: list of (x,y) or 'inf'."""
    return evaluate(plane_to_sphere(pts), **kw)


def best_planar_picture(ev, k=0):
    """Planar coordinates (stereographic projection from the k-th best centre)."""
    deg, o = ev.best_centres(k + 1)[k]
    return sphere_to_plane(ev.V, o), deg, o


if __name__ == "__main__":
    # sanity checks on the known records
    tests = {
        "n=6 record (8)": [(0, 0), (20, 0), (5, 15), (5, 0), (2, 6), (10, 10)],
        "n=7 record (11)": [(0, 0), (20, 0), (5, 15), (5, 0), (2, 6), (10, 10), (5, 5)],
        "n=8 record (17)": [(0, 0), (0, 5), (0, 10), (0, 15), (3, 6), (3, 9), (5, 5), (5, 10)],
        "two squares (18)": [(1, 0), (0, 1), (-1, 0), (0, -1), (2, 0), (0, 2), (-2, 0), (0, -2)],
        "orthocentric+feet (n=7) + inf": [(0, 3), (1, 0), (3, 0), (0, -1), (1.2, -0.6), (2, 1), (0, 0), 'inf'],
    }
    for name, pts in tests.items():
        ev = evaluate_planar(pts)
        s = ev.summary()
        # planar count with O = north pole (= Euclidean count when 'inf' not in pts)
        north = np.array([0, 0, 1.0])
        eu = ev.nb - len(blocks_through(north, ev.Nn, ev.D))
        print(name, s, "euclid", eu)
