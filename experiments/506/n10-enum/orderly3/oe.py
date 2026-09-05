#!/usr/bin/env python3
"""
oe.py -- orderly enumeration of abstract Moebius block structures (Erdos #506), written from
scratch (orderly3 agent, 2026-09-05).  Nothing is imported from the earlier enumerators.

Setting.  P = {0..n-1}.  A *rich block* is a subset of size 4..n-1 (a circle or line through >= 4
points).  F = family of rich blocks, pairwise sharing <= 2 points.  A *line set* L is a family of
lines (blocks through infinity): rich blocks of F used as lines, or *uncovered* triples (triples
contained in no rich block), pairwise sharing <= 1 point.
   circles = C(n,3) - D(F) - |L|,   D(F) = sum_{B in F} (C(|B|,3) - 1).
Caps (necessary conditions from Sylvester-Gallai on the derived point sets):
   (C2) for every p: sum_{B in F, p in B} C(|B|-1, 2) <= capD = C(n-1,2) - o(n-1)
        and #{4-blocks through p} <= max4
   (C3) sum_{L in lines} C(|L|,2) <= capL = C(n,2) - o(n); #3-lines <= t3n; |lines| <= ellmax.
Mode 'sg' uses o(m) = 1 only (Sylvester-Gallai); mode 'table' uses the known o(m) and orchard
numbers t3(m) (cited: Kelly-Moser / Csima-Sawyer / Crowe-McKee for o, Burr-Gruenbaum-Sloane 1974
for t3).

Canonical form of a family: the lexicographically smallest sorted tuple of block bitmasks over
all relabellings (S_n).  Two-phase orderly generation:
  phase 1: canonical families of *big* blocks (size >= 5), full S_n canonicity via the "marked
           first block" test (equivalent to brute force, see is_canonical_marked);
  phase 2: for each canonical F_big, orderly extension by 4-blocks, canonicity of the 4-block list
           tested against Aut(F_big) (canonical form of the pair = big list first, then 4-list).
Parent-of-canonical-is-canonical (removing the largest block) holds for both phases, so the DFS
(adding blocks in increasing mask order, keeping only canonical children) visits every
isomorphism class exactly once.  Validated against brute force in validate.py.
"""
import sys, json, time, itertools, math
from math import comb
import numpy as np

# ----------------------------------------------------------------------------------------------
# tables
O_TABLE = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5, 11: 6, 12: 6, 13: 6, 14: 7}   # min #ordinary lines
T3_TABLE = {3: 1, 4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10, 10: 12, 11: 16, 12: 19}               # orchard numbers
INF = 10 ** 9


def popcount(m):
    return bin(m).count('1')


def mask_to_list(m, n):
    return [i for i in range(n) if m >> i & 1]


def list_to_mask(s):
    return sum(1 << i for i in s)


class Caps:
    def __init__(self, n, mode):
        self.n, self.mode = n, mode
        if mode == 'none':
            self.capD, self.max4, self.capL, self.t3n, self.ellmax = INF, INF, INF, INF, INF
            return
        if mode == 'sg':
            o_prev, o_n, t3_prev, t3_n = 1, 1, INF, INF
        elif mode == 'table':
            o_prev, o_n, t3_prev, t3_n = O_TABLE[n - 1], O_TABLE[n], T3_TABLE[n - 1], T3_TABLE[n]
        else:
            raise ValueError(mode)
        self.capD = comb(n - 1, 2) - o_prev
        self.max4 = min(self.capD // 3, t3_prev)
        self.capL = comb(n, 2) - o_n
        self.t3n = t3_n
        # ellmax: l lines cover >= 3*min(l,t3n) + 6*max(0, l - t3n) pairs
        l = 0
        while 3 * min(l + 1, self.t3n) + 6 * max(0, l + 1 - self.t3n) <= self.capL:
            l += 1
        self.ellmax = l

    def __str__(self):
        return (f"Caps(n={self.n}, mode={self.mode}, capD={self.capD}, max4={self.max4}, "
                f"capL={self.capL}, t3n={self.t3n}, ellmax={self.ellmax})")


# ----------------------------------------------------------------------------------------------
# permutations and images
def all_perms(n):
    return np.array(list(itertools.permutations(range(n))), dtype=np.int8)


_stab_cache = {}


def stab_perms(n, k):
    """All permutations of {0..n-1} fixing {0..k-1} setwise (S_k x S_{n-k}), as an (N, n) array
    with row sigma meaning i -> sigma[i]."""
    key = (n, k)
    if key not in _stab_cache:
        a = list(itertools.permutations(range(k)))
        b = list(itertools.permutations(range(k, n)))
        _stab_cache[key] = np.array([p + q for p in a for q in b], dtype=np.int8)
    return _stab_cache[key]


def images(perms, masks):
    """perms: (N, n) int8 (row = permutation i -> perms[r, i]); masks: list of ints.
    Returns (N, len(masks)) int32: image mask of each block under each permutation."""
    n = perms.shape[1]
    pw = (np.int32(1) << perms.astype(np.int32))          # (N, n): 2^{sigma(i)}
    out = np.empty((perms.shape[0], len(masks)), dtype=np.int32)
    for j, m in enumerate(masks):
        idx = [i for i in range(n) if m >> i & 1]
        out[:, j] = pw[:, idx].sum(axis=1)
    return out


def rows_less_eq(rows, target):
    """rows: (N, k) already sorted per row; target: length-k sequence.
    Returns (some_row_is_lex_smaller, boolean mask of rows equal to target)."""
    N = rows.shape[0]
    lt = np.zeros(N, dtype=bool)
    eq = np.ones(N, dtype=bool)
    for j, t in enumerate(target):
        col = rows[:, j]
        lt |= eq & (col < t)
        eq &= (col == t)
    return bool(lt.any()), eq


def apply_perm_mask(perm, m):
    """perm: sequence with perm[i] = image of i."""
    r = 0
    i = 0
    while m:
        if m & 1:
            r |= 1 << perm[i]
        m >>= 1
        i += 1
    return r


def canon_bruteforce(n, fam, perms=None):
    """Lexicographically smallest sorted mask tuple over all of S_n (validation only)."""
    if perms is None:
        perms = all_perms(n)
    if len(fam) == 0:
        return ()
    rows = np.sort(images(perms, list(fam)), axis=1)
    order = np.lexsort(rows.T[::-1])            # primary key = column 0
    return tuple(int(x) for x in rows[order[0]])


def is_canonical_marked(n, fam):
    """Full-S_n canonicity test of a sorted family of masks (sizes may differ), plus Aut(fam).
    A permutation sigma with sorted(sigma(fam)) <= fam must send some block of minimum size k to
    {0..k-1} (the smallest possible mask of a block of size >= k), so sigma = tau o rho_c with
    rho_c a fixed permutation sending block c to {0..k-1} and tau in the setwise stabiliser of
    {0..k-1}.  Every such sigma is examined exactly once.
    Returns (is_canonical, aut) with aut an (A, n) int8 array of all automorphisms
    (aut is None when not canonical)."""
    fam = tuple(fam)
    if len(fam) == 0:
        return True, all_perms(n)
    k = min(popcount(m) for m in fam)
    if fam[0] != (1 << k) - 1:
        return False, None
    stab = stab_perms(n, k)
    auts = []
    for c in fam:
        if popcount(c) != k:
            continue
        pts = mask_to_list(c, n)
        rest = [i for i in range(n) if not (c >> i & 1)]
        rho = [0] * n
        for j, p in enumerate(pts):
            rho[p] = j
        for j, p in enumerate(rest):
            rho[p] = k + j
        fam2 = sorted(apply_perm_mask(rho, m) for m in fam)
        rows = np.sort(images(stab, fam2), axis=1)
        less, eq = rows_less_eq(rows, fam)
        if less:
            return False, None
        taus = stab[eq]                       # tau with sorted(tau(rho(fam))) == fam
        auts.append(taus[:, rho])             # sigma(i) = tau(rho(i))
    return True, np.concatenate(auts, axis=0)


# ----------------------------------------------------------------------------------------------
# phase 1: canonical families of big blocks
def phase1(n, sizes, caps, verbose=True):
    """Orderly generation of all canonical families of blocks with sizes in `sizes` (all >= 5),
    pairwise sharing <= 2 points, satisfying cap (C2).  Returns list of (fam_tuple, aut_array)."""
    big_masks = []
    for s in sizes:
        big_masks += [list_to_mask(c) for c in itertools.combinations(range(n), s)]
    big_masks.sort()
    out = []
    tests = [0]

    def used_of(fam):
        u = [0] * n
        for B in fam:
            w = comb(popcount(B) - 1, 2)
            for p in mask_to_list(B, n):
                u[p] += w
        return u

    def rec(fam, aut):
        out.append((fam, aut))
        last = fam[-1] if fam else -1
        used = used_of(fam)
        for b in big_masks:
            if b <= last:
                continue
            if any(popcount(b & B) > 2 for B in fam):
                continue
            w = comb(popcount(b) - 1, 2)
            if any(used[p] + w > caps.capD for p in mask_to_list(b, n)):
                continue
            child = tuple(sorted(fam + (b,)))
            tests[0] += 1
            ok, a = is_canonical_marked(n, child)
            if ok:
                rec(child, a)

    rec((), all_perms(n) if n <= 8 else None)
    if verbose:
        print(f"phase 1: {len(out)} canonical big-block families (incl. empty), "
              f"{tests[0]} canonicity tests", flush=True)
    return out


# ----------------------------------------------------------------------------------------------
# line sets
def uncovered_triples(n, F):
    tr = []
    for c in itertools.combinations(range(n), 3):
        m = list_to_mask(c)
        if not any(m & B == m for B in F):
            tr.append(m)
    return tr


def line_sets(n, F, caps, need, want_all=True, limit=200000):
    """All sets L of lines (rich blocks of F or uncovered triples), pairwise sharing <= 1 point,
    with sum C(|L|,2) <= capL, #3-lines <= t3n, |L| >= need.  Returns (lmax, list_of_sets)."""
    cands = [(m, popcount(m)) for m in F] + [(m, 3) for m in uncovered_triples(n, F)]
    cands.sort(key=lambda t: (-t[1], t[0]))          # rich first
    K = len(cands)
    masks = [c[0] for c in cands]
    sizes = [c[1] for c in cands]
    compat = [[popcount(masks[i] & masks[j]) <= 1 for j in range(K)] for i in range(K)]
    best = [0]
    found = []

    def rec(start, chosen, pairs, n3, avail):
        L = len(chosen)
        if L > best[0]:
            best[0] = L
        if L >= need and need > 0:
            found.append(tuple(chosen))
            if len(found) > limit:
                raise RuntimeError("too many line sets")
        # bound
        if L + len(avail) < need:
            return
        for ii, i in enumerate(avail):
            s = sizes[i]
            if pairs + comb(s, 2) > caps.capL:
                continue
            if s == 3 and n3 + 1 > caps.t3n:
                continue
            new_avail = [j for j in avail[ii + 1:] if compat[i][j]]
            if L + 1 + len(new_avail) < need:
                continue
            chosen.append(masks[i])
            rec(i + 1, chosen, pairs + comb(s, 2), n3 + (s == 3), new_avail)
            chosen.pop()

    rec(0, [], 0, 0, list(range(K)))
    if need <= 0:
        found = [()]
    return best[0], found


# ----------------------------------------------------------------------------------------------
# phase 2: extension by 4-blocks
class Phase2:
    """Orderly extension of a canonical big-block family F_big by 4-blocks.

    Codes of 4-blocks = their index in the lexicographically sorted list of 4-tuples.
    Canonical form of (F_big, F4) = (sorted big masks, sorted 4-codes) compared big part first;
    with F_big canonical this holds iff sorted(sigma(F4)) >= F4 for all sigma in Aut(F_big).
    Only automorphisms sending some block of F4 to a code <= F4[0] can violate the inequality
    ("relevant" automorphisms; the set grows monotonically along a DFS path and is maintained
    incrementally).  Upper bounds on the number of 4-blocks that can still be added (each block
    uses 3 pairs at each of its 4 points and 4 uncovered triples):
      a_p = min( (capD - used_p)//3, max4 - d4_p, #candidates containing p,
                 floor( sum_q floor(thirds_pq / 2) / 3 ) )
      where thirds_pq = #{x : {p,q,x} is a triple of some remaining candidate} (blocks through
      {p,q} have pairwise disjoint remaining pairs);  #additional <= floor(sum_p a_p / 4).
    Candidates are used in increasing code order, so the bound for the j-th child only counts the
    candidates from position j on (suffix bound): the loop breaks when it drops below need."""

    def __init__(self, n, Fbig, aut, caps, need_total, log=None, collect_only=False):
        self.n, self.Fbig, self.caps, self.need_total = n, tuple(Fbig), caps, need_total
        self.collect_only = collect_only
        self.collected = []
        # codes of 4-blocks: lexicographic order of the sorted point tuples (all blocks through
        # point 0 first, then those through 1 but not 0, ...): closes the points one by one along
        # the DFS, which makes the per-point saturation bounds effective.
        self.masks4 = [list_to_mask(c) for c in itertools.combinations(range(n), 4)]
        self.M = len(self.masks4)
        self.pts4 = np.array([mask_to_list(m, n) for m in self.masks4], dtype=np.int8)   # (M,4)
        self.onehot = np.zeros((self.M, n), dtype=np.int16)
        for j in range(self.M):
            self.onehot[j, self.pts4[j]] = 1
        self.mask_arr = np.array(self.masks4, dtype=np.int32)
        self.index_of = np.full(1 << n, -1, dtype=np.int32)
        self.index_of[self.mask_arr] = np.arange(self.M, dtype=np.int32)
        # triples and pairs
        self.triples = [list_to_mask(c) for c in itertools.combinations(range(n), 3)]
        tri_index = {m: i for i, m in enumerate(self.triples)}
        self.tri_of_block = np.array([[tri_index[m & ~(1 << p)] for p in mask_to_list(m, n)]
                                      for m in self.masks4], dtype=np.int32)          # (M,4)
        pairs = list(itertools.combinations(range(n), 2))
        self.pair_tri = np.array([[tri_index[(1 << a) | (1 << b) | (1 << x)] for x in range(n) if x != a and x != b]
                                  for (a, b) in pairs], dtype=np.int32)                # (P, n-2)
        self.pt_pair = np.zeros((n, len(pairs)), dtype=np.int64)
        for k, (a, b) in enumerate(pairs):
            self.pt_pair[a, k] = 1
            self.pt_pair[b, k] = 1
        # automorphism images as codes (drop the identity)
        if aut is None:
            raise ValueError("phase 2 needs Aut(F_big)")
        ident = np.all(aut == np.arange(n, dtype=np.int8), axis=1)
        aut_nt = aut[~ident]
        self.A = aut_nt.shape[0]
        self.idx_img = self.index_of[images(aut_nt, self.masks4)] if self.A else None   # (A, M)
        # initial per-point usage and compatibility with F_big
        self.used0 = np.zeros(n, dtype=np.int32)
        self.D0 = 0
        for B in self.Fbig:
            s = popcount(B)
            self.D0 += comb(s, 3) - 1
            for p in mask_to_list(B, n):
                self.used0[p] += comb(s - 1, 2)
        pc = np.array([popcount(m) for m in range(1 << n)], dtype=np.int8)
        ok = np.ones(self.M, dtype=bool)
        for B in self.Fbig:
            ok &= pc[self.mask_arr & B] <= 2
        self.compat0 = ok
        self.pc = pc
        self.Dmin = need_total - caps.ellmax
        self.nodes = self.evals = self.tests = 0
        self.results = []
        self.log = log

    # ----- bounds
    def point_caps(self, cands, used, d4):
        """a_p as in the class docstring, for the candidate list cands (array of codes)."""
        caps = self.caps
        a = np.minimum((caps.capD - used) // 3, caps.max4 - d4).astype(np.int64)
        if len(cands) == 0:
            return np.zeros(self.n, dtype=np.int64)
        tri_av = np.zeros(len(self.triples), dtype=np.int64)
        tri_av[self.tri_of_block[cands].ravel()] = 1
        thirds = tri_av[self.pair_tri].sum(axis=1)                  # (P,)
        per_pt = self.pt_pair @ (thirds // 2)                        # (n,)
        a = np.minimum(a, per_pt // 3)
        a = np.minimum(a, self.onehot[cands].sum(axis=0))
        return np.maximum(a, 0)

    def bound_extra(self, cands, used, d4):
        if len(cands) == 0:
            return 0
        return min(len(cands), int(self.point_caps(cands, used, d4).sum()) // 4)

    # ----- canonicity
    def canonical(self, F4, rel):
        """rel: boolean array over the non-identity automorphisms (relevant ones)."""
        if self.A == 0:
            return True
        idx = np.flatnonzero(rel)
        if len(idx) == 0:
            return True
        self.tests += 1
        rows = np.sort(self.idx_img[np.ix_(idx, F4)], axis=1)
        less, _ = rows_less_eq(rows, F4)
        return not less

    def evaluate(self, F4):
        self.evals += 1
        F = list(self.Fbig) + [self.masks4[i] for i in F4]
        D = self.D0 + 3 * len(F4)
        need = self.need_total - D
        lmax, sets = line_sets(self.n, F, self.caps, need)
        if lmax >= need:
            cnt = comb(self.n, 3) - D - lmax
            deg = [sum(1 for B in F if B >> p & 1) for p in range(self.n)]
            rec = {"n": self.n, "mode": self.caps.mode,
                   "blocks": [mask_to_list(B, self.n) for B in F],
                   "sizes": [popcount(B) for B in F], "D": D, "ell_max": lmax, "count_min": cnt,
                   "degrees": deg,
                   "line_sets": [[mask_to_list(L, self.n) for L in S] for S in sets]}
            self.results.append(rec)
            if self.log:
                print(f"    CANDIDATE: sizes={rec['sizes']} D={D} ell_max={lmax} count={cnt} "
                      f"degrees={deg} #line_sets={len(sets)}", file=self.log, flush=True)

    def run(self):
        caps = self.caps
        n = self.n
        used = self.used0.copy()
        d4 = np.zeros(n, dtype=np.int32)
        point_ok = (used + 3 <= caps.capD) & (d4 < caps.max4)
        feas = self.compat0 & point_ok[self.pts4].all(axis=1)
        cands = np.nonzero(feas)[0]
        self.root_bound = self.bound_extra(cands, used, d4)
        if not self.collect_only and self.D0 + 3 * self.root_bound + caps.ellmax < self.need_total:
            return self.results
        rel = np.zeros(self.A, dtype=bool)
        self.t0 = time.time()
        self._rec([], cands, used, d4, rel)
        return self.results

    def _rec(self, F4, cands, used, d4, rel):
        caps = self.caps
        self.nodes += 1
        if self.log and self.nodes % 50000 == 0:
            print(f"    ... nodes={self.nodes} tests={self.tests} evals={self.evals} depth={len(F4)} "
                  f"[{time.time()-self.t0:.0f}s]", file=self.log, flush=True)
        D = self.D0 + 3 * len(F4)
        if self.collect_only:
            self.collected.append(tuple(sorted(list(self.Fbig) + [self.masks4[i] for i in F4])))
            need_extra = 0
        else:
            if D >= self.Dmin:
                self.evaluate(F4)
            need_extra = max(0, -(-(self.need_total - caps.ellmax - D) // 3))   # ceil
        if len(cands) == 0:
            return
        pts4, pc, mask_arr = self.pts4, self.pc, self.mask_arr
        if need_extra > 0:
            a = self.point_caps(cands, used, d4)
            suffix = np.cumsum(self.onehot[cands][::-1], axis=0)[::-1]        # (m, n)
            sbound = np.minimum(suffix, a).sum(axis=1) // 4                   # (m,)
        c0 = F4[0] if F4 else None
        for pos in range(len(cands)):
            if need_extra > 0 and sbound[pos] < need_extra:
                break
            c = int(cands[pos])
            used2 = used.copy()
            d42 = d4.copy()
            for p in pts4[c]:
                used2[p] += 3
                d42[p] += 1
            point_ok = (used2 + 3 <= caps.capD) & (d42 < caps.max4)
            rest = cands[pos + 1:]
            keep = (pc[mask_arr[rest] & mask_arr[c]] <= 2) & point_ok[pts4[rest]].all(axis=1)
            rest = rest[keep]
            if need_extra > 1 and self.bound_extra(rest, used2, d42) < need_extra - 1:
                continue
            child = F4 + [c]
            first = c0 if c0 is not None else c
            rel2 = (rel | (self.idx_img[:, c] <= first)) if self.A else rel
            if not self.canonical(child, rel2):
                continue
            self._rec(child, rest, used2, d42, rel2)


# ----------------------------------------------------------------------------------------------
def run_enum(n, target, mode, sizes, outfile=None, skip_empty=True, verbose=True):
    caps = Caps(n, mode)
    need_total = comb(n, 3) - target
    print(f"n={n} target={target} mode={mode} need D+ell>={need_total} {caps} sizes={sizes}",
          flush=True)
    t0 = time.time()
    fams = phase1(n, sizes, caps)
    print(f"phase 1 done in {time.time()-t0:.1f}s", flush=True)
    all_results = []
    total_nodes = 0
    for fi, (fam, aut) in enumerate(fams):
        if len(fam) == 0 and skip_empty:
            print(f"[{fi}] empty big family: skipped (all-4-block case handled separately)", flush=True)
            continue
        if aut is None:
            ok, aut = is_canonical_marked(n, fam)
            assert ok
        t1 = time.time()
        ph = Phase2(n, fam, aut, caps, need_total, log=sys.stdout if verbose else None)
        res = ph.run()
        total_nodes += ph.nodes
        print(f"[{fi}] sizes={[popcount(B) for B in fam]} blocks={[mask_to_list(B, n) for B in fam]} "
              f"D0={ph.D0} |Aut|={ph.A+1} root_bound={ph.root_bound} nodes={ph.nodes} "
              f"evals={ph.evals} candidates={len(res)} [{time.time()-t1:.1f}s]", flush=True)
        all_results += res
    print(f"TOTAL phase-2 nodes={total_nodes}; candidate structures (count <= {target}): "
          f"{len(all_results)} [{time.time()-t0:.1f}s]", flush=True)
    for i, r in enumerate(all_results):
        print(f"  cand {i}: sizes={r['sizes']} D={r['D']} ell_max={r['ell_max']} count={r['count_min']} "
              f"degrees={r['degrees']}\n     blocks={r['blocks']}")
    if outfile:
        with open(outfile, 'w') as f:
            json.dump({"n": n, "target": target, "mode": mode, "caps": str(caps), "sizes": sizes,
                       "results": all_results}, f, indent=1)
    return all_results


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('n', type=int)
    ap.add_argument('target', type=int)
    ap.add_argument('--mode', default='sg')
    ap.add_argument('--sizes', default=None, help='comma list of big block sizes (default 5..n-1)')
    ap.add_argument('--out', default=None)
    ap.add_argument('--with-empty', action='store_true', help='also run the all-4-block case '
                    '(phase 2 with Aut = S_n; only sensible for n <= 9)')
    a = ap.parse_args()
    sizes = [int(s) for s in a.sizes.split(',')] if a.sizes else list(range(5, a.n))
    run_enum(a.n, a.target, a.mode, sizes, outfile=a.out, skip_empty=not a.with_empty)
