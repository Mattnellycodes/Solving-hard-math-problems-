"""Local-search engine for Erdős #506 (search-local agent, 2026-09-05).

Universe U = finite planar points (floats; exact data kept separately) plus optionally the point at
infinity (always the LAST index when present).  Blocks of U = circles through >= 3 finite points, and
lines through >= 2 finite points (a line is the block through infinity; it is a block of a subset S
only when |line ∩ S| >= 3, counting infinity).  Every triple of U lies in exactly one block, so for
S ⊆ U the Möbius blocks of S are {B : |B ∩ S| >= 3} and

    best_count(S) = |B(S)| - max_{O ∉ S} deg_S(O),

where deg_S(O) = number of blocks of S through O.  In the search O ranges over universe points
(including infinity) and optional extra "O-only" points; the true optimum over all O (pairwise
intersections of blocks of S) is computed by `true_count` for the final candidates.
"""
import math, itertools, time, sys
import numpy as np
from numba import njit

BIG = 10 ** 6


# ---------------------------------------------------------------------------- block finding
def _quantise(vals, q):
    return np.round(vals / q).astype(np.int64)


def find_blocks_float(P, q=1e-7, has_inf=True, chunk=2_000_000):
    """P: (N,2) float array of finite points.  Returns list of blocks (sorted index tuples) and a
    parallel list of geometric descriptors ('C', cx, cy, r) / ('L', A, B, C) (line: A x + B y = C,
    (A,B) unit).  Lines through >= 2 finite points are included (with index N = infinity appended)
    if has_inf, otherwise only lines through >= 3 finite points (no infinity index)."""
    N = len(P)
    scale = max(1.0, float(np.abs(P).max()))
    x = P[:, 0] / scale; y = P[:, 1] / scale          # normalise to [-1,1]
    circ = {}    # key -> set of points
    circ_geom = {}
    lines = {}
    line_geom = {}
    # lines from pairs (exact enough: normal vector quantised)
    I, J = np.triu_indices(N, 1)
    A = y[J] - y[I]; B = x[I] - x[J]
    nrm = np.hypot(A, B)
    A = A / nrm; B = B / nrm; C = A * x[I] + B * y[I]
    flip = (A < -1e-12) | ((np.abs(A) <= 1e-12) & (B < 0))
    A[flip] *= -1; B[flip] *= -1; C[flip] *= -1
    kA = _quantise(A, q); kB = _quantise(B, q); kC = _quantise(C, q)
    for t in range(len(I)):
        key = (kA[t], kB[t], kC[t])
        s = lines.get(key)
        if s is None:
            lines[key] = {int(I[t]), int(J[t])}; line_geom[key] = (A[t] * 1.0, B[t] * 1.0, C[t] * 1.0)
        else:
            s.add(int(I[t])); s.add(int(J[t]))
    lines = _merge_neighbours(lines, line_geom)
    # circles from triples
    ntrip = N * (N - 1) * (N - 2) // 6
    trip_iter = itertools.combinations(range(N), 3)
    done = 0
    while done < ntrip:
        T = np.fromiter(itertools.chain.from_iterable(itertools.islice(trip_iter, chunk)), dtype=np.int64)
        if len(T) == 0:
            break
        T = T.reshape(-1, 3)
        done += len(T)
        ax, ay = x[T[:, 0]], y[T[:, 0]]; bx, by = x[T[:, 1]], y[T[:, 1]]; cx, cy = x[T[:, 2]], y[T[:, 2]]
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        # collinear if the triangle area is tiny relative to the side lengths
        side = np.maximum(np.hypot(ax - bx, ay - by), np.maximum(np.hypot(ax - cx, ay - cy), np.hypot(bx - cx, by - cy)))
        col = np.abs(d) < 1e-9 * side * side
        ok = ~col
        a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx * cx + cy * cy
        dd = np.where(ok, d, 1.0)
        ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / dd
        uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / dd
        r = np.hypot(ax - ux, ay - uy)
        # circles with huge radius relative to the universe are treated as circles anyway (the
        # collinearity test above is the only line test); quantise centre and radius
        ku = _quantise(ux, q); kv = _quantise(uy, q); kr = _quantise(r, q)
        idx = np.nonzero(ok)[0]
        for t in idx:
            key = (ku[t], kv[t], kr[t])
            s = circ.get(key)
            if s is None:
                circ[key] = {int(T[t, 0]), int(T[t, 1]), int(T[t, 2])}
                circ_geom[key] = (ux[t] * 1.0, uy[t] * 1.0, r[t] * 1.0)
            else:
                s.add(int(T[t, 0])); s.add(int(T[t, 1])); s.add(int(T[t, 2]))
    circ = _merge_neighbours(circ, circ_geom)
    blocks, geom = [], []
    for key, s in circ.items():
        g = circ_geom[key]
        blocks.append(tuple(sorted(s))); geom.append(('C', g[0] * scale, g[1] * scale, g[2] * scale))
    for key, s in lines.items():
        g = line_geom[key]
        if has_inf:
            blocks.append(tuple(sorted(s)) + (N,)); geom.append(('L', g[0], g[1], g[2] * scale))
        elif len(s) >= 3:
            blocks.append(tuple(sorted(s))); geom.append(('L', g[0], g[1], g[2] * scale))
    return blocks, geom


def find_blocks_exact(pts, has_inf=True):
    """Exact block finder for rational points (list of (Fraction, Fraction)); Python integers.
    Returns (blocks, geom) in the same format as find_blocks_float (geometry as floats)."""
    from fractions import Fraction as Fr
    from math import gcd
    N = len(pts)
    den = 1
    for x, y in pts:
        den = den * x.denominator // gcd(den, x.denominator)
        den = den * y.denominator // gcd(den, y.denominator)
    X = [int(x * den) for x, y in pts]; Y = [int(y * den) for x, y in pts]
    assert len(set(zip(X, Y))) == N, "duplicate points"
    S2 = [X[i] * X[i] + Y[i] * Y[i] for i in range(N)]
    circ = {}; lines = {}
    for i in range(N):
        xi, yi, si = X[i], Y[i], S2[i]
        for j in range(i + 1, N):
            xj, yj = X[j], Y[j]
            dx1 = xj - xi; dy1 = yj - yi; s1 = S2[j] - si
            # line key for the pair (used for 2-point lines when has_inf)
            A = dy1; B = -dx1; C = A * xi + B * yi
            g = gcd(gcd(abs(A), abs(B)), abs(C)) or 1
            A //= g; B //= g; C //= g
            if A < 0 or (A == 0 and B < 0):
                A, B, C = -A, -B, -C
            lk = (A, B, C)
            s = lines.get(lk)
            if s is None:
                lines[lk] = {i, j}
            else:
                s.add(i); s.add(j)
            for k in range(j + 1, N):
                xk, yk = X[k], Y[k]
                dx2 = xk - xi; dy2 = yk - yi
                det = dx1 * dy2 - dx2 * dy1
                if det == 0:
                    continue   # collinear: handled by the pair lines
                s2 = S2[k] - si
                un = s1 * dy2 - s2 * dy1; vn = dx1 * s2 - dx2 * s1; d = 2 * det
                g = gcd(gcd(abs(un), abs(vn)), abs(d))
                un //= g; vn //= g; d //= g
                if d < 0:
                    un, vn, d = -un, -vn, -d
                r2n = (xi * d - un) ** 2 + (yi * d - vn) ** 2
                key = (un, vn, d, r2n)
                s = circ.get(key)
                if s is None:
                    circ[key] = {i, j, k}
                else:
                    s.add(i); s.add(j); s.add(k)
    blocks, geom = [], []
    for (un, vn, d, r2n), s in circ.items():
        blocks.append(tuple(sorted(s)))
        geom.append(('C', un / d / den, vn / d / den, math.sqrt(r2n) / d / den))
    for (A, B, C), s in lines.items():
        nrm = math.hypot(A, B)
        if has_inf:
            blocks.append(tuple(sorted(s)) + (N,)); geom.append(('L', A / nrm, B / nrm, C / nrm / den))
        elif len(s) >= 3:
            blocks.append(tuple(sorted(s))); geom.append(('L', A / nrm, B / nrm, C / nrm / den))
    return blocks, geom


def _merge_neighbours(groups, geom):
    """merge groups whose quantised keys differ by at most 1 in every coordinate (boundary effects)."""
    keys = list(groups.keys())
    parent = {k: k for k in keys}

    def find(k):
        while parent[k] != k:
            parent[k] = parent[parent[k]]; k = parent[k]
        return k
    keyset = set(keys)
    for k in keys:
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    if dx == 0 and dy == 0 and dz == 0:
                        continue
                    k2 = (k[0] + dx, k[1] + dy, k[2] + dz)
                    if k2 in keyset:
                        a, b = find(k), find(k2)
                        if a != b:
                            parent[a] = b
    out = {}
    for k in keys:
        r = find(k)
        if r not in out:
            out[r] = set()
        out[r] |= groups[k]
    return out


def check_blocks(P, blocks, geom, has_inf, tol=1e-6):
    """sanity: every member lies on its block; every triple in exactly one block (count check)."""
    N = len(P)
    scale = max(1.0, float(np.abs(P).max()))
    for b, g in zip(blocks, geom):
        fin = [i for i in b if i < N]
        if g[0] == 'C':
            _, cx, cy, r = g
            for i in fin:
                assert abs(math.hypot(P[i, 0] - cx, P[i, 1] - cy) - r) < tol * scale, (b, g)
        else:
            _, A, B, C = g
            for i in fin:
                assert abs(A * P[i, 0] + B * P[i, 1] - C) < tol * scale, (b, g)
    M = N + (1 if has_inf else 0)
    tot = sum(len(b) * (len(b) - 1) * (len(b) - 2) // 6 for b in blocks)
    assert tot == M * (M - 1) * (M - 2) // 6, (tot, M)
    return True


# ---------------------------------------------------------------------------- CSR structure
class Universe:
    def __init__(self, name, P, has_inf=True, exact=None, labels=None, q=1e-7, extra_O=None):
        """P: (N,2) float array of finite points.  exact: list of exact coordinate pairs (Fractions or
        sympy) parallel to P (None if not available).  extra_O: (K,2) float array of O-only points."""
        self.name = name
        self.P = np.asarray(P, dtype=np.float64)
        self.N = len(self.P)
        self.has_inf = has_inf
        self.exact = exact
        self.labels = labels
        self.M = self.N + (1 if has_inf else 0)          # selectable points
        t = time.time()
        from fractions import Fraction as _Fr
        if exact is not None and len(exact) == self.N and all(isinstance(e[0], _Fr) and isinstance(e[1], _Fr) for e in exact):
            blocks, geom = find_blocks_exact(exact, has_inf=has_inf)
            self.block_mode = 'exact'
        else:
            blocks, geom = find_blocks_float(self.P, q=q, has_inf=has_inf)
            self.block_mode = 'float'
        check_blocks(self.P, blocks, geom, has_inf, tol=1e-6 if self.block_mode == 'float' else 1e-3)
        self.blocks = blocks; self.geom = geom
        self.t_build = time.time() - t
        # O-only extra points: append after the selectable ones; they belong to blocks they lie on
        self.K = 0
        if extra_O is not None and len(extra_O):
            self.add_extra_O(np.asarray(extra_O, dtype=np.float64))
        self._csr()

    def add_extra_O(self, E, tol=1e-7):
        scale = max(1.0, float(np.abs(self.P).max()))
        self.extra = E
        self.K = len(E)
        newblocks = []
        for b, g in zip(self.blocks, self.geom):
            if g[0] == 'C':
                on = np.abs(np.hypot(E[:, 0] - g[1], E[:, 1] - g[2]) - g[3]) < tol * scale
            else:
                on = np.abs(g[1] * E[:, 0] + g[2] * E[:, 1] - g[3]) < tol * scale
            idx = np.nonzero(on)[0]
            newblocks.append(tuple(b) + tuple(int(self.M + i) for i in idx))
        self.blocks = newblocks

    def _csr(self):
        T = self.M + self.K
        nb = len(self.blocks)
        bptr = np.zeros(nb + 1, np.int64)
        for i, b in enumerate(self.blocks):
            bptr[i + 1] = bptr[i] + len(b)
        bmem = np.zeros(bptr[-1], np.int64)
        for i, b in enumerate(self.blocks):
            bmem[bptr[i]:bptr[i + 1]] = b
        deg = np.zeros(T, np.int64)
        for b in self.blocks:
            for p in b:
                deg[p] += 1
        pptr = np.zeros(T + 1, np.int64); pptr[1:] = np.cumsum(deg)
        pblk = np.zeros(pptr[-1], np.int64)
        fill = pptr[:-1].copy()
        for i, b in enumerate(self.blocks):
            for p in b:
                pblk[fill[p]] = i; fill[p] += 1
        self.bptr, self.bmem, self.pptr, self.pblk = bptr, bmem, pptr, pblk
        self.T = T
        # sizes (selectable points only)
        self.sizes = {}
        for b in self.blocks:
            k = sum(1 for p in b if p < self.M)
            self.sizes[k] = self.sizes.get(k, 0) + 1

    def label(self, i):
        if i == self.N and self.has_inf:
            return 'inf'
        if i >= self.M:
            return f'O{i - self.M}'
        if self.labels is not None:
            return self.labels[i]
        return f'({self.P[i, 0]:.6g},{self.P[i, 1]:.6g})'

    def summary(self):
        return (f"{self.name}: N={self.N} inf={self.has_inf} extraO={self.K} blocks={len(self.blocks)} mode={getattr(self, 'block_mode', '?')} "
                f"sizes={dict(sorted(self.sizes.items()))} build={self.t_build:.1f}s")


# ---------------------------------------------------------------------------- numba kernels
@njit(cache=True)
def k_add(p, n, inS, cnt, deg, bptr, bmem, pptr, pblk):
    inS[p] = True
    dnb = 0; dfull = 0
    for t in range(pptr[p], pptr[p + 1]):
        b = pblk[t]
        cnt[b] += 1
        c = cnt[b]
        if c == 3:
            dnb += 1
            for u in range(bptr[b], bptr[b + 1]):
                deg[bmem[u]] += 1
        if c == n:
            dfull += 1
    return dnb, dfull


@njit(cache=True)
def k_remove(p, n, inS, cnt, deg, bptr, bmem, pptr, pblk):
    inS[p] = False
    dnb = 0; dfull = 0
    for t in range(pptr[p], pptr[p + 1]):
        b = pblk[t]
        c = cnt[b]
        cnt[b] = c - 1
        if c == 3:
            dnb -= 1
            for u in range(bptr[b], bptr[b + 1]):
                deg[bmem[u]] -= 1
        if c == n:
            dfull -= 1
    return dnb, dfull


@njit(cache=True)
def k_value(nb, full, inS, deg, T):
    """objective = nb - max deg over points not in S (all T points incl. O-only); BIG if degenerate"""
    if full > 0:
        return BIG, -1
    md = -1; mo = -1
    for p in range(T):
        if not inS[p] and deg[p] > md:
            md = deg[p]; mo = p
    return nb - md, mo


@njit(cache=True)
def k_eval(S, n, T, bptr, bmem, pptr, pblk):
    inS = np.zeros(T, np.bool_)
    cnt = np.zeros(bptr.shape[0] - 1, np.int64)
    deg = np.zeros(T, np.int64)
    nb = 0; full = 0
    for i in range(n):
        a, b = k_add(S[i], n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
    v, o = k_value(nb, full, inS, deg, T)
    return v, o, nb


@njit(cache=True)
def k_sa(n, M, T, bptr, bmem, pptr, pblk, iters, T0, T1, seed, init, tabu_len, bias):
    """Simulated annealing over n-subsets of the M selectable points (0..M-1); O ranges over all T
    points.  Moves: swap one point out, one in.  Candidate to add is chosen with bias towards high
    current degree (bias = number of random candidates sampled, best-degree one taken) half the time.
    Returns best value, best subset, best centre."""
    np.random.seed(seed)
    inS = np.zeros(T, np.bool_)
    cnt = np.zeros(bptr.shape[0] - 1, np.int64)
    deg = np.zeros(T, np.int64)
    S = init.copy()
    nb = 0; full = 0
    for i in range(n):
        a, b = k_add(S[i], n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
    cur, co = k_value(nb, full, inS, deg, T)
    best = cur; bestS = S.copy(); bestO = co
    tabu = np.full(max(tabu_len, 1), -1, np.int64); tpos = 0
    lam = math.log(T1 / T0) / max(iters - 1, 1)
    for it in range(iters):
        temp = T0 * math.exp(lam * it)
        # pick point to remove
        ri = np.random.randint(n)
        p_out = S[ri]
        # pick point to add
        if bias > 1 and np.random.random() < 0.5:
            bestd = -1; p_in = -1
            for _ in range(bias):
                q = np.random.randint(M)
                if inS[q]:
                    continue
                tb = False
                for z in range(tabu.shape[0]):
                    if tabu[z] == q:
                        tb = True
                if tb:
                    continue
                if deg[q] > bestd:
                    bestd = deg[q]; p_in = q
            if p_in < 0:
                continue
        else:
            p_in = np.random.randint(M)
            if inS[p_in]:
                continue
            tb = False
            for z in range(tabu.shape[0]):
                if tabu[z] == p_in:
                    tb = True
            if tb:
                continue
        a, b = k_remove(p_out, n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
        a, b = k_add(p_in, n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
        new, no = k_value(nb, full, inS, deg, T)
        delta = new - cur
        if delta <= 0 or np.random.random() < math.exp(-delta / temp):
            S[ri] = p_in
            cur = new; co = no
            if tabu_len > 0:
                tabu[tpos] = p_out; tpos = (tpos + 1) % tabu.shape[0]
            if cur < best:
                best = cur; bestS[:] = S; bestO = co
        else:
            a, b = k_remove(p_in, n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
            a, b = k_add(p_out, n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
    return best, bestS, bestO


@njit(cache=True)
def k_exhaustive(n, M, T, bptr, bmem, pptr, pblk, thresh, einf):
    """Enumerate all n-subsets of the M selectable points (lexicographic DFS).  Returns
    (total, best value, histogram of values (<400), list of subsets with value <= thresh (max 20000))."""
    inS = np.zeros(T, np.bool_)
    cnt = np.zeros(bptr.shape[0] - 1, np.int64)
    deg = np.zeros(T, np.int64)
    S = np.empty(n, np.int64)
    hist = np.zeros(400, np.int64)
    hist_e = np.zeros(400, np.int64)
    hist_nb = np.zeros(400, np.int64)
    good = np.zeros((20000, n), np.int64); ngood = 0
    best = BIG; total = 0
    nb = 0; full = 0
    depth = 0; S[0] = -1
    while depth >= 0:
        S[depth] += 1
        if S[depth] > M - (n - depth):
            depth -= 1
            if depth >= 0:
                a, b = k_remove(S[depth], n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
            continue
        a, b = k_add(S[depth], n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
        if depth == n - 1:
            total += 1
            if full == 0:
                v, o = k_value(nb, full, inS, deg, T)
                if v < 400:
                    hist[v] += 1
                if einf >= 0:
                    e = nb - deg[einf]   # Euclidean count (O = infinity)
                    if e < 400:
                        hist_e[e] += 1
                if nb < 400:
                    hist_nb[nb] += 1
                if v < best:
                    best = v
                if v <= thresh and ngood < 20000:
                    good[ngood, :] = S; ngood += 1
            a, b = k_remove(S[depth], n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
        else:
            depth += 1
            S[depth] = S[depth - 1]
    return total, best, hist, hist_e, hist_nb, good[:ngood]


# ---------------------------------------------------------------------------- true objective
def blocks_of_subset(U, S):
    """Möbius blocks of S (as index tuples restricted to S) with their geometry."""
    Sset = set(int(i) for i in S)
    out = []
    for b, g in zip(U.blocks, U.geom):
        m = [i for i in b if i in Sset]
        if len(m) >= 3:
            out.append((tuple(m), g))
    return out


def _intersections(g1, g2):
    """intersection points (floats) of two blocks given by geometry descriptors."""
    out = []
    if g1[0] == 'L' and g2[0] == 'L':
        _, A1, B1, C1 = g1; _, A2, B2, C2 = g2
        det = A1 * B2 - A2 * B1
        if abs(det) < 1e-12:
            return out
        out.append(((C1 * B2 - C2 * B1) / det, (A1 * C2 - A2 * C1) / det))
        return out
    if g1[0] == 'L':
        g1, g2 = g2, g1
    _, cx, cy, r = g1
    if g2[0] == 'L':
        _, A, B, C = g2
        d = A * cx + B * cy - C
        if abs(d) > r * (1 + 1e-9) + 1e-12:
            return out
        h2 = r * r - d * d
        h = math.sqrt(h2) if h2 > 0 else 0.0
        px, py = cx - A * d, cy - B * d
        out.append((px - B * h, py + A * h)); out.append((px + B * h, py - A * h))
        return out
    _, dx, dy, r2 = g2
    D = math.hypot(dx - cx, dy - cy)
    if D < 1e-12 or D > (r + r2) * (1 + 1e-9) + 1e-12 or D < abs(r - r2) * (1 - 1e-9) - 1e-12:
        return out
    a = (r * r - r2 * r2 + D * D) / (2 * D); h2 = r * r - a * a
    h = math.sqrt(h2) if h2 > 0 else 0.0
    mx, my = cx + a * (dx - cx) / D, cy + a * (dy - cy) / D
    out.append((mx + h * (dy - cy) / D, my - h * (dx - cx) / D)); out.append((mx - h * (dy - cy) / D, my + h * (dx - cx) / D))
    return out


def true_count(U, S, tol=1e-9, excl=1e-4):
    """|B(S)| - max_O deg_S(O) over ALL O (pairwise intersections of blocks of S, plus infinity,
    plus universe/extra points), O not in S.  Returns (count, O, deg, nblocks, lines).
    Residual tolerance tol*scale (intersections are double precision, so genuine concurrences have
    residuals ~1e-12*scale); candidates within excl*scale of a point of S are rejected (a centre that
    close to a point of S is a numerical artifact: circles through that point pass near it).
    A value below the combinatorial proxy must be certified exactly before it is trusted."""
    S = [int(i) for i in S]
    Sset = set(S)
    BS = blocks_of_subset(U, S)
    nb = len(BS)
    scale = max(1.0, float(np.abs(U.P).max()))
    finite_pts = [U.P[i] for i in S if i < U.N]
    inf_in_S = U.has_inf and (U.N in Sset)
    lines = sum(1 for m, g in BS if g[0] == 'L')
    best_deg = 0 if inf_in_S else lines
    best_O = None if inf_in_S else 'inf'

    def deg_at(x, y):
        d = 0
        for m, g in BS:
            if g[0] == 'C':
                if abs(math.hypot(x - g[1], y - g[2]) - g[3]) < tol * scale:
                    d += 1
            else:
                if abs(g[1] * x + g[2] * y - g[3]) < tol * scale:
                    d += 1
        return d

    def near_S(x, y):
        return any(abs(x - p[0]) < excl * scale and abs(y - p[1]) < excl * scale for p in finite_pts)
    cands = {}
    for (m1, g1), (m2, g2) in itertools.combinations(BS, 2):
        for (x, y) in _intersections(g1, g2):
            if near_S(x, y):
                continue
            key = (round(x / (1e-5 * scale)), round(y / (1e-5 * scale)))
            cands[key] = (x, y)
    for (x, y) in cands.values():
        d = deg_at(x, y)
        if d > best_deg:
            best_deg = d; best_O = (x, y)
    return nb - best_deg, best_O, best_deg, nb, lines


def formula(n):
    return (n - 1) * (n - 2) // 2 + 1 - (n - 1) // 2


def describe(U, S):
    BS = blocks_of_subset(U, S)
    sizes = {}
    for m, g in BS:
        sizes[len(m)] = sizes.get(len(m), 0) + 1
    return dict(sorted(sizes.items()))


# ---------------------------------------------------------------------------- triple-table kernels
# For universes with M selectable points, tab[(i*M+j)*M+k] (i<j<k) = index of the unique block containing
# the triple.  Adding/removing a point p costs O(n^2) table lookups instead of O(#blocks through p).
def build_triple_table(U):
    M = U.M
    tab = np.full(M * M * M, -1, dtype=np.int32)
    _fill_triple_table(tab, M, U.bptr, U.bmem)
    assert (tab[np.array([(i * M + j) * M + k for i in range(M) for j in range(i + 1, M) for k in range(j + 1, M)])] >= 0).all()
    return tab


@njit(cache=True)
def _fill_triple_table(tab, M, bptr, bmem):
    nb = bptr.shape[0] - 1
    for b in range(nb):
        lo, hi = bptr[b], bptr[b + 1]
        for u in range(lo, hi):
            i = bmem[u]
            if i >= M:
                continue
            for v in range(u + 1, hi):
                j = bmem[v]
                if j >= M:
                    continue
                for w in range(v + 1, hi):
                    k = bmem[w]
                    if k >= M:
                        continue
                    a, bb, c = i, j, k
                    if a > bb:
                        a, bb = bb, a
                    if bb > c:
                        bb, c = c, bb
                    if a > bb:
                        a, bb = bb, a
                    tab[(a * M + bb) * M + c] = b


@njit(cache=True)
def _blocks_via(p, S, n, M, tab, stamp, stampval, blist, bmult):
    """distinct blocks through p and >= 2 other points of S (S given as array of n entries; p may or
    may not be in S).  Returns count; blist[:cnt] = block ids, bmult[:cnt] = number of pairs found
    (= C(k,2) where k = number of other S-points on the block)."""
    cnt = 0
    for a in range(n):
        q = S[a]
        if q == p:
            continue
        for b2 in range(a + 1, n):
            r = S[b2]
            if r == p:
                continue
            i, j, k = p, q, r
            if i > j:
                i, j = j, i
            if j > k:
                j, k = k, j
            if i > j:
                i, j = j, i
            b = tab[(i * M + j) * M + k]
            if stamp[b] != stampval:
                stamp[b] = stampval
                blist[cnt] = b; bmult[cnt] = 1; cnt += 1
            else:
                for t in range(cnt):
                    if blist[t] == b:
                        bmult[t] += 1
                        break
    return cnt


@njit(cache=True)
def _k_from_pairs(m):
    # m = C(k,2) -> k
    k = 2
    while k * (k - 1) // 2 < m:
        k += 1
    return k


@njit(cache=True)
def t_add(p, S, m, n, M, tab, stamp, stampval, blist, bmult, inS, deg, bptr, bmem):
    """add p to S (S[:m] holds the m current points, p not among them); n = target subset size for
    the degeneracy check; updates deg; returns (dnb, dfull)"""
    inS[p] = True
    cnt = _blocks_via(p, S, m, M, tab, stamp, stampval, blist, bmult)
    dnb = 0; dfull = 0
    for t in range(cnt):
        b = blist[t]
        k = _k_from_pairs(bmult[t]) + 1      # new number of S-points on b
        if k == 3:
            dnb += 1
            for u in range(bptr[b], bptr[b + 1]):
                deg[bmem[u]] += 1
        if k == n:
            dfull += 1
    return dnb, dfull


@njit(cache=True)
def t_remove(p, S, m, n, M, tab, stamp, stampval, blist, bmult, inS, deg, bptr, bmem):
    """remove p from S (S[:m] holds the m current points incl. p); returns (dnb, dfull)"""
    inS[p] = False
    cnt = _blocks_via(p, S, m, M, tab, stamp, stampval, blist, bmult)
    dnb = 0; dfull = 0
    for t in range(cnt):
        b = blist[t]
        k = _k_from_pairs(bmult[t]) + 1      # old number of S-points on b
        if k == 3:
            dnb -= 1
            for u in range(bptr[b], bptr[b + 1]):
                deg[bmem[u]] -= 1
        if k == n:
            dfull -= 1
    return dnb, dfull


@njit(cache=True)
def t_eval(S, n, M, T, tab, bptr, bmem, nblocks):
    inS = np.zeros(T, np.bool_)
    deg = np.zeros(T, np.int64)
    stamp = np.zeros(nblocks, np.int64)
    blist = np.zeros(n * n, np.int64); bmult = np.zeros(n * n, np.int64)
    cur = np.zeros(n, np.int64)
    nb = 0; full = 0; sv = 0
    for i in range(n):
        cur[i] = S[i]
        sv += 1
        a, b = t_add(S[i], cur, i, n, M, tab, stamp, sv, blist, bmult, inS, deg, bptr, bmem); nb += a; full += b
    v, o = k_value(nb, full, inS, deg, T)
    return v, o, nb


@njit(cache=True)
def t_sa(n, M, T, tab, bptr, bmem, nblocks, iters, T0, T1, seed, init, tabu_len, bias):
    """SA with the triple table (same moves/acceptance as k_sa)."""
    np.random.seed(seed)
    inS = np.zeros(T, np.bool_)
    deg = np.zeros(T, np.int64)
    stamp = np.zeros(nblocks, np.int64)
    blist = np.zeros(n * n, np.int64); bmult = np.zeros(n * n, np.int64)
    S = init.copy()
    nb = 0; full = 0; sv = 0
    for i in range(n):
        sv += 1
        a, b = t_add(S[i], S, i, n, M, tab, stamp, sv, blist, bmult, inS, deg, bptr, bmem); nb += a; full += b
    cur, co = k_value(nb, full, inS, deg, T)
    best = cur; bestS = S.copy(); bestO = co
    tabu = np.full(max(tabu_len, 1), -1, np.int64); tpos = 0
    lam = math.log(T1 / T0) / max(iters - 1, 1)
    for it in range(iters):
        temp = T0 * math.exp(lam * it)
        ri = np.random.randint(n)
        p_out = S[ri]
        if bias > 1 and np.random.random() < 0.5:
            bestd = -1; p_in = -1
            for _ in range(bias):
                q = np.random.randint(M)
                if inS[q]:
                    continue
                tb = False
                for z in range(tabu.shape[0]):
                    if tabu[z] == q:
                        tb = True
                if tb:
                    continue
                if deg[q] > bestd:
                    bestd = deg[q]; p_in = q
            if p_in < 0:
                continue
        else:
            p_in = np.random.randint(M)
            if inS[p_in]:
                continue
            tb = False
            for z in range(tabu.shape[0]):
                if tabu[z] == p_in:
                    tb = True
            if tb:
                continue
        # remove p_out: S currently holds n points incl. p_out
        sv += 1
        a, b = t_remove(p_out, S, n, n, M, tab, stamp, sv, blist, bmult, inS, deg, bptr, bmem); nb += a; full += b
        # move p_out to the end and add p_in in its place logically: S[:n-1] must be the others
        S[ri] = S[n - 1]; S[n - 1] = p_out
        sv += 1
        a, b = t_add(p_in, S, n - 1, n, M, tab, stamp, sv, blist, bmult, inS, deg, bptr, bmem); nb += a; full += b
        S[n - 1] = p_in
        new, no = k_value(nb, full, inS, deg, T)
        delta = new - cur
        if delta <= 0 or np.random.random() < math.exp(-delta / temp):
            cur = new; co = no
            if tabu_len > 0:
                tabu[tpos] = p_out; tpos = (tpos + 1) % tabu.shape[0]
            if cur < best:
                best = cur; bestS[:] = S; bestO = co
        else:
            sv += 1
            a, b = t_remove(p_in, S, n, n, M, tab, stamp, sv, blist, bmult, inS, deg, bptr, bmem); nb += a; full += b
            sv += 1
            a, b = t_add(p_out, S, n - 1, n, M, tab, stamp, sv, blist, bmult, inS, deg, bptr, bmem); nb += a; full += b
            S[n - 1] = p_out
    return best, bestS, bestO
