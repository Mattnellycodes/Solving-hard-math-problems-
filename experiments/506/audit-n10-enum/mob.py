"""
Independent audit enumerator for Erdos #506 abstract Moebius block structures (audit-n10-enum).
Written from scratch; nothing imported from other agents' code.

Model (n points 0..n-1):
  rich blocks = subsets of size 4..n-1, pairwise sharing <= 2 points;
  D = sum over rich blocks of (C(k,3)-1);
  lines L = subset of (rich blocks  U  triples not covered by any rich block), pairwise sharing <= 1 point,
            sum over lines of C(|l|,2) <= capL   (Sylvester-Gallai at infinity: capL = C(n,2)-o(n));
  per point p: sum over rich blocks B through p of C(|B|-1,2) <= capD  (derived Sylvester-Gallai: capD = C(n-1,2)-o(n-1));
  circles = C(n,3) - D - |L|.
Goal: all isomorphism classes of (F, L) with circles <= target, i.e. D + |L| >= need = C(n,3) - target.

Isomorph rejection: stage 1 (blocks of size >= 5) canonical form = lexicographically least sorted bitmask list
under all of S_n (brute force over n! permutations with numpy); stage 2 (4-blocks) canonical = lexicographically
least sorted list of 4-tuple codes under Aut(F_big) (returned by the brute force).  Both are Read-Faradzev orderly
generations (parent = remove largest code); validated against brute force in validate.py.
"""
import numpy as np, itertools, math, sys, time
from math import comb

def all_perms(n):
    P = np.zeros((1, 1), dtype=np.int8)
    for k in range(2, n + 1):
        m = P.shape[0]
        out = np.empty((m * k, k), dtype=np.int8)
        for pos in range(k):
            blk = out[pos * m:(pos + 1) * m]
            blk[:, :pos] = P[:, :pos]
            blk[:, pos] = k - 1
            blk[:, pos + 1:] = P[:, pos:]
        P = out
    return P

def popcount(x):
    return bin(x).count('1')

def bits(m):
    out = []
    i = 0
    while m:
        if m & 1:
            out.append(i)
        m >>= 1
        i += 1
    return out

def mask_of(s):
    r = 0
    for i in s:
        r |= 1 << i
    return r

# ---------------------------------------------------------------- stage 1 machinery
class PermTable:
    """All permutations of range(n) with their 'shift' table 1<<perm[i]."""
    def __init__(self, n):
        self.n = n
        self.P = all_perms(n)
        self.N = self.P.shape[0]
        self.shift = (np.int16(1) << self.P.astype(np.int16))  # (N, n)

    def image_masks(self, masks):
        """(N, k) int16 array: sorted image masks of the family under every permutation."""
        k = len(masks)
        imgs = np.empty((self.N, k), dtype=np.int16)
        for j, m in enumerate(masks):
            pts = bits(m)
            imgs[:, j] = self.shift[:, pts].sum(axis=1, dtype=np.int16)
        imgs.sort(axis=1)
        return imgs

def rows_less_or_equal(imgs, target):
    """Return (any_row_less, equal_rows_boolean) comparing each row lexicographically with target."""
    N = imgs.shape[0]
    diff = imgs.astype(np.int32) - np.asarray(target, dtype=np.int32)[None, :]
    nz = diff != 0
    has = nz.any(axis=1)
    first = nz.argmax(axis=1)
    val = diff[np.arange(N), first]
    less = np.any(has & (val < 0))
    equal = ~has
    return less, equal

def stage1(n, capD, sizes, PT=None, verbose=True):
    """All canonical families of blocks of sizes in `sizes` (>=5), pairwise sharing <= 2 points, derived cap capD.
    Returns list of (masks sorted, Aut as int8 array (|Aut|, n))."""
    if PT is None:
        PT = PermTable(n)
    bigs = sorted(mask_of(c) for k in sizes for c in itertools.combinations(range(n), k))
    results = []
    tests = [0]
    def canonical(F):
        # necessary: smallest mask must be {0..k-1} for the smallest size present
        kmin = min(popcount(m) for m in F)
        if F[0] != (1 << kmin) - 1:
            return None
        tests[0] += 1
        imgs = PT.image_masks(F)
        less, equal = rows_less_or_equal(imgs, F)
        if less:
            return None
        return PT.P[equal]
    def rec(F, cov, aut):
        results.append((list(F), aut))
        last = F[-1] if F else -1
        for m in bigs:
            if m <= last:
                continue
            if any(popcount(m & b) > 2 for b in F):
                continue
            k = popcount(m)
            add = comb(k - 1, 2)
            if any(cov[p] + add > capD for p in bits(m)):
                continue
            F2 = sorted(F + [m])
            aut2 = canonical(F2)
            if aut2 is None:
                continue
            cov2 = cov[:]
            for p in bits(m):
                cov2[p] += add
            rec(F2, cov2, aut2)
    rec([], [0] * n, PT.P)
    if verbose:
        print(f"stage1: n={n} capD={capD} sizes={list(sizes)}: {len(results)} canonical families (incl. empty), {tests[0]} brute-force tests", flush=True)
    return results

# ---------------------------------------------------------------- lines
def line_sets(n, rich_masks, capL, need_lines, limit=10**6):
    """All sets of lines (subsets of rich blocks U uncovered triples, pairwise sharing <= 1 point,
    total pairs <= capL) with at least need_lines lines.  Returns (ell_max, list of line sets as masks)."""
    covered = set()
    for m in rich_masks:
        for t in itertools.combinations(bits(m), 3):
            covered.add(mask_of(t))
    pool = list(rich_masks) + [mask_of(t) for t in itertools.combinations(range(n), 3) if mask_of(t) not in covered]
    w = [comb(popcount(m), 2) for m in pool]
    K = len(pool)
    compat = [0] * K
    for i in range(K):
        for j in range(K):
            if i != j and popcount(pool[i] & pool[j]) <= 1:
                compat[i] |= 1 << j
    best = [0]
    out = []
    def rec(cur, avail, weight):
        if len(cur) > best[0]:
            best[0] = len(cur)
        if len(cur) >= need_lines:
            out.append([pool[i] for i in cur])
            if len(out) > limit:
                raise RuntimeError("line set limit")
        if len(cur) + popcount(avail) < need_lines:
            return
        a = avail
        while a:
            i = (a & -a).bit_length() - 1
            a &= a - 1
            if weight + w[i] > capL:
                continue
            # remaining after i: only indices > i
            rest = avail & compat[i] & ~((1 << (i + 1)) - 1)
            if len(cur) + 1 + popcount(rest) < need_lines:
                # later i have even fewer successors? not monotone in general, so just continue
                continue
            rec(cur + [i], rest, weight + w[i])
    rec([], (1 << K) - 1, 0)
    return best[0], out

# ---------------------------------------------------------------- stage 2
def code4(t):
    return t[0] * 1000 + t[1] * 100 + t[2] * 10 + t[3]

class Stage2:
    def __init__(self, n, fixed_masks, G, capD, capL, need, degcap=None, ellmax=None, verbose=False, record_all_lines=True):
        """fixed_masks: rich blocks already present (any sizes); G: int8 array of permutations (automorphisms of the
        fixed family) under which 4-block families are canonicalised; need = required D + ell."""
        self.n = n
        self.fixed = list(fixed_masks)
        self.G = np.ascontiguousarray(G)
        self.capD, self.capL, self.need = capD, capL, need
        self.ellmax = capL // 3 if ellmax is None else ellmax
        self.degcap = [capD // 3] * n if degcap is None else list(degcap)
        self.verbose = verbose
        self.D0 = sum(comb(popcount(m), 3) - 1 for m in self.fixed)
        cov = [0] * n
        deg = [0] * n
        for m in self.fixed:
            k = popcount(m)
            for p in bits(m):
                cov[p] += comb(k - 1, 2)
                if k == 4:
                    deg[p] += 1
        self.cov0, self.deg0 = cov, deg
        cands = []
        for c in itertools.combinations(range(n), 4):
            m = mask_of(c)
            if any(popcount(m & b) > 2 for b in self.fixed):
                continue
            if any(cov[p] + 3 > capD or deg[p] + 1 > self.degcap[p] for p in c):
                continue
            cands.append(c)
        cands.sort()
        self.cands = cands
        self.cmask = [mask_of(c) for c in cands]
        self.ccode = [code4(c) for c in cands]
        K = len(cands)
        self.K = K
        self.compat = [0] * K
        for i in range(K):
            for j in range(K):
                if i != j and popcount(self.cmask[i] & self.cmask[j]) <= 2:
                    self.compat[i] |= 1 << j
        self.pmask = [0] * n
        for i, c in enumerate(cands):
            for p in c:
                self.pmask[p] |= 1 << i
        self.nodes = 0
        self.evals = 0
        self.candidates = []
        self.record_all_lines = record_all_lines
        self.collect_only = False

    def bound(self, avail, cov, deg):
        """Upper bound on the number of further 4-blocks that can be added from the candidate set `avail`."""
        n = self.n
        T = {}
        cnt = [0] * n
        a = avail
        while a:
            i = (a & -a).bit_length() - 1
            a &= a - 1
            c = self.cands[i]
            m = self.cmask[i]
            for p in c:
                cnt[p] += 1
                for q in c:
                    if q > p:
                        T[(p, q)] = T.get((p, q), 0) | (m & ~(1 << p) & ~(1 << q))
        thirds = [0] * n
        for (p, q), tm in T.items():
            h = popcount(tm) // 2
            thirds[p] += h
            thirds[q] += h
        tot = 0
        for p in range(n):
            ap = min(cnt[p], (self.capD - cov[p]) // 3, self.degcap[p] - deg[p], thirds[p] // 3)
            tot += max(ap, 0)
        return tot // 4

    def canonical(self, F4idx):
        if self.G.shape[0] <= 1 or not F4idx:
            return True
        pts = np.array([self.cands[i] for i in F4idx], dtype=np.int64)  # (m,4)
        codes = np.array([self.ccode[i] for i in F4idx], dtype=np.int32)
        img = self.G[:, pts].astype(np.int32)  # (N, m, 4)
        img.sort(axis=2)
        ic = img[..., 0] * 1000 + img[..., 1] * 100 + img[..., 2] * 10 + img[..., 3]
        ic.sort(axis=1)
        less, _ = rows_less_or_equal(ic, codes)
        return not less

    def evaluate(self, F4idx, D):
        self.evals += 1
        rich = self.fixed + [self.cmask[i] for i in F4idx]
        need_lines = max(self.need - D, 0)
        ellmax, sets = line_sets(self.n, rich, self.capL, need_lines)
        if D + ellmax >= self.need:
            degs = [0] * self.n
            for m in rich:
                for p in bits(m):
                    degs[p] += 1
            rec = dict(blocks=[bits(m) for m in rich], D=D, ell_max=ellmax,
                       count_min=comb(self.n, 3) - D - ellmax, degrees=degs,
                       line_sets=[[bits(m) for m in s] for s in sets] if self.record_all_lines else None)
            self.candidates.append(rec)
            if self.verbose:
                print(f"    CANDIDATE D={D} ell_max={ellmax} count={rec['count_min']} degrees={degs} #line_sets={len(sets)}", flush=True)

    def run(self):
        avail0 = (1 << self.K) - 1
        self.rec([], avail0, self.cov0[:], self.deg0[:], self.D0)
        return self.candidates

    def rec(self, F4idx, avail, cov, deg, D):
        self.nodes += 1
        if self.collect_only:
            self.candidates.append(sorted(self.fixed + [self.cmask[i] for i in F4idx]))
        elif D + self.ellmax >= self.need:
            self.evaluate(F4idx, D)
        a = avail
        while a:
            i = (a & -a).bit_length() - 1
            a &= a - 1
            suffix = a  # candidates with index > i
            if D + 3 * (1 + self.bound(suffix, cov, deg)) + self.ellmax < self.need:
                break  # monotone in i
            c = self.cands[i]
            cov2 = cov[:]
            deg2 = deg[:]
            child = suffix & self.compat[i]
            for p in c:
                cov2[p] += 3
                deg2[p] += 1
                if cov2[p] + 3 > self.capD or deg2[p] + 1 > self.degcap[p]:
                    child &= ~self.pmask[p]
            if D + 3 + 3 * self.bound(child, cov2, deg2) + self.ellmax < self.need:
                continue
            F2 = F4idx + [i]
            if not self.canonical(F2):
                continue
            self.rec(F2, child, cov2, deg2, D + 3)

# ---------------------------------------------------------------- full pipeline (big-block case)
def caps(n, mode):
    if mode == 'sg':
        capD = comb(n - 1, 2) - 1
        capL = comb(n, 2) - 1
    else:
        o = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5, 11: 6, 12: 6, 13: 6, 14: 7}
        capD = comb(n - 1, 2) - o[n - 1]
        capL = comb(n, 2) - o[n]
    return capD, capL

def run_big(n, target, mode='sg', sizes=None, verbose=True, PT=None):
    capD, capL = caps(n, mode)
    need = comb(n, 3) - target
    if sizes is None:
        sizes = range(5, n)
    t0 = time.time()
    fams = stage1(n, capD, sizes, PT=PT, verbose=verbose)
    if verbose:
        print(f"stage1 time {time.time()-t0:.1f}s", flush=True)
    allc = []
    tot_nodes = 0
    for idx, (F, aut) in enumerate(fams):
        if not F:
            continue
        t1 = time.time()
        S = Stage2(n, F, aut, capD, capL, need, verbose=verbose)
        cands = S.run()
        tot_nodes += S.nodes
        if verbose:
            print(f"[{idx}] big={[bits(m) for m in F]} |Aut|={aut.shape[0]} D0={S.D0} #cands4={S.K} nodes={S.nodes} evals={S.evals} candidates={len(cands)} [{time.time()-t1:.1f}s]", flush=True)
        for c in cands:
            c['big'] = [bits(m) for m in F]
        allc.extend(cands)
    if verbose:
        print(f"TOTAL stage-2 nodes={tot_nodes}; candidate structures (count <= {target}): {len(allc)} [{time.time()-t0:.1f}s]", flush=True)
    return allc

# ---------------------------------------------------------------- canonical forms of full structures (for dedupe / comparison)
def canon_full(PT, blocks):
    """Brute-force S_n canonical form of a family of blocks (any sizes): lexicographically least sorted mask list."""
    masks = sorted(mask_of(b) for b in blocks)
    imgs = PT.image_masks(masks)
    # lexicographic minimum row
    order = np.lexsort(imgs.T[::-1])
    return tuple(int(x) for x in imgs[order[0]])

def canon_full_with_lines(PT, blocks, lines):
    """Canonical form of (F, L): lines are marked by adding bit n (so a rich block used as a line is a different
    element from the same block not used as a line; 3-lines are separate elements)."""
    n = PT.n
    masks = sorted([mask_of(b) for b in blocks] + [mask_of(l) | (1 << n) for l in lines])
    k = len(masks)
    imgs = np.empty((PT.N, k), dtype=np.int32)
    for j, m in enumerate(masks):
        pts = bits(m & ((1 << n) - 1))
        imgs[:, j] = PT.shift[:, pts].sum(axis=1, dtype=np.int32) + (m & (1 << n))
    imgs.sort(axis=1)
    order = np.lexsort(imgs.T[::-1])
    return tuple(int(x) for x in imgs[order[0]])
