#!/usr/bin/env python3
"""orderly2.py -- two-phase orderly (Read-Faradzev) generation of abstract Moebius block structures
for Erdos #506.  Written from scratch by the n10-enum agent (2026-09-05); independent of
../theory/mobius_enum*.py and of the CP-SAT models.

Objects.  P = {0..n-1}.  A structure is a family F of rich blocks (subsets of size 4..n-1) such that two
rich blocks share <= 2 points (C1).  Every triple lies in exactly one block (the triples in no rich block
are the 3-blocks).  Deficit D(F) = sum_{B in F} (C(|B|,3) - 1);  circles(P) = C(n,3) - D(F) - l, where
l = number of lines = blocks through infinity: a set of blocks (rich blocks of F or 3-blocks) pairwise
sharing <= 1 point.

Necessary conditions used (theorems; see ../theory/REPORT.md Lemma 2.2/2.8 and ../verify_independent):
 (C2)  derived Sylvester-Gallai cap: for every p, sum_{B in F, p in B} C(|B|-1,2) <= C(n-1,2) - o(n-1),
       because inversion about p maps the blocks through p to the >=3-point lines of a non-collinear real
       (n-1)-point set, which has >= o(n-1) ordinary (2-point) lines.
 (C2') orchard cap: the number of 4-blocks through p (= derived lines with exactly 3 points) <= t3(n-1).
 (C3)  line caps: sum_L C(|L|,2) <= C(n,2) - o(n) and the number of 3-point lines <= t3(n).
 o(m)  = minimum number of ordinary lines of m non-collinear points: o(m) >= 1 (Sylvester-Gallai); exact
       values o(3..14) = 3,3,4,3,3,4,6,5,6,6,6,7 (Crowe-McKee 1968; OEIS A003034; Pach-Sharir table).
 t3(m) = orchard number = max number of 3-point lines of m points: t3(3..12) = 1,1,2,4,6,7,10,12,16,19
       (Burr-Gruenbaum-Sloane 1974; OEIS A003035).  Our use only needs "lines with exactly 3 points are
       at most t3", which holds under either reading of the definition.
Modes: 'table' = exact o-table + t3;  'sg_t3' = o = 1 + t3;  'sg' = o = 1 only;  'none' = no caps
(used for validation against brute force).

Canonical forms (my own design, validated in validate.py):
 Phase 1: families F_big of blocks of size >= 5; code of a block = its bitmask; F_big is canonical iff the
   sorted code list is lexicographically least over S_n.  Test: the least code of a canonical family is
   2^k - 1 (k = smallest block size), so only permutations mapping a smallest block onto {0..k-1} can
   produce a smaller or equal sorted list; there are k!(n-k)! of them per smallest block and we test them
   all (vectorised).  Aut(F_big) is read off the same computation.
 Phase 2: 4-blocks with code = rank in lexicographic order of sorted 4-tuples (so all blocks through
   point 0 come first, which makes the per-point saturation prune effective).  (F_big, F4) is canonical
   iff F_big is canonical (phase 1) and sorted(sigma(F4)) >= sorted(F4) for all sigma in Aut(F_big).
   Only automorphisms mapping some block of F4 onto the least code of F4 can violate this, and this set
   is maintained incrementally (rel).
 Orderly generation: children append a code larger than all present codes.  If a child is canonical so
   is its parent (lex-min sorted lists: removing the largest element of a lex-least list leaves a
   lex-least list -- if sigma made the parent smaller, inserting sigma(b) could only make the child
   smaller as well).  Hence a DFS over canonical families with children in increasing code order visits
   every isomorphism class of pairs (F_big, F4) exactly once.

Pruning in phase 2 (all upper bounds on the number of 4-blocks that can still be added):
   per point p: rem_p = min((cap_derived - cov_p)//3, t3 - deg4_p, (sum_q floor(thirds_pq/2))//3), where
   thirds_pq = number of x with {p,q,x} uncovered and contained in a remaining candidate (blocks through
   {p,q} pairwise share only {p,q}, so their third/fourth points are disjoint pairs);
   availability: candidates come in increasing code order, so only the suffix R[j:] is available, and
   #additional <= sum_p min(rem_p, #{c in R[j:] : p in c}) // 4.
   Prune when D + 3 * #additional < D_min = C(n,3) - target - ell_max.

Usage: python3 orderly2.py n target [--mode table|sg_t3|sg|none] [--big 5,6,...] [--include-empty]
                                     [--families i,j,...] [--out results.json]
"""
import sys, itertools, math, time, json, argparse
from math import comb
import numpy as np

O_TABLE = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5, 11: 6, 12: 6, 13: 6, 14: 7}
T3 = {3: 1, 4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10, 10: 12, 11: 16, 12: 19}
INF = 10 ** 9


def popcount(x):
    return bin(x).count('1')


def ceil_div(a, b):
    return -((-a) // b)


class Params:
    def __init__(self, n, target, mode='table', big_sizes=None):
        assert mode in ('table', 'sg_t3', 'sg', 'none')
        self.n, self.target, self.mode = n, target, mode
        if mode == 'none':
            self.cap_derived = INF; self.cap_lines = INF; self.max4 = INF; self.max3lines = INF
            self.ell_max = INF
        else:
            o = (lambda m: O_TABLE[m]) if mode == 'table' else (lambda m: 1)
            self.cap_derived = comb(n - 1, 2) - o(n - 1)
            self.cap_lines = comb(n, 2) - o(n)
            use_t3 = mode != 'sg'
            self.max4 = min(self.cap_derived // 3, T3[n - 1] if use_t3 else INF)
            self.max3lines = T3[n] if use_t3 else INF
            t = min(self.max3lines, self.cap_lines // 3)
            self.ell_max = t + (self.cap_lines - 3 * t) // 6      # 3-lines are cheapest, then 4-lines
        self.N3 = comb(n, 3)
        self.D_min = self.N3 - target - self.ell_max if self.ell_max < INF else -INF
        self.big_sizes = list(big_sizes) if big_sizes is not None else list(range(5, n))

    def describe(self):
        return (f"n={self.n} target={self.target} mode={self.mode} C(n,3)={self.N3} cap_derived={self.cap_derived} "
                f"max4/pt={self.max4} cap_lines={self.cap_lines} max3lines={self.max3lines} ell_max={self.ell_max} "
                f"D_min={self.D_min} big_sizes={self.big_sizes}")


def rows_less(M, ref):
    """Boolean array: which rows of M (2-D, each row sorted, same length as ref) are lexicographically < ref."""
    neq = M != ref
    has = neq.any(axis=1)
    first = neq.argmax(axis=1)
    vals = M[np.arange(M.shape[0]), first]
    return has & (vals < ref[first])


# ----------------------------------------------------------------------------------------------------
# Phase 1: canonical families of big blocks (size >= 5) under S_n
# ----------------------------------------------------------------------------------------------------
class Phase1:
    def __init__(self, P):
        self.P = P
        n = self.n = P.n
        self.size = [popcount(m) for m in range(1 << n)]
        self.pts = [[i for i in range(n) if m >> i & 1] for m in range(1 << n)]
        self.big = [m for m in range(1 << n) if self.size[m] in P.big_sizes]     # increasing code order
        self.perm_cache = {}
        self.families = []
        self.nodes = 0

    def block_perms(self, mask):
        """All permutations sigma (rows of an (m, n) int8 array, sigma[i] = new label of point i)
        with sigma(block) = {0..k-1}, k = |block|."""
        k = self.size[mask]
        if k not in self.perm_cache:
            A = np.array(list(itertools.permutations(range(k))), dtype=np.int8).reshape(-1, k)
            C = np.array(list(itertools.permutations(range(k, self.n))), dtype=np.int8).reshape(-1, self.n - k)
            self.perm_cache[k] = (A, C)
        A, C = self.perm_cache[k]
        B = self.pts[mask]
        R = [i for i in range(self.n) if not mask >> i & 1]
        S = np.empty((len(A) * len(C), self.n), dtype=np.int8)
        S[:, B] = np.repeat(A, len(C), axis=0)
        if R:
            S[:, R] = np.tile(C, (len(A), 1))
        return S

    def relevant_perms(self, F):
        k = min(self.size[m] for m in F)
        return np.concatenate([self.block_perms(m) for m in F if self.size[m] == k], axis=0)

    def images(self, F, S):
        """(m, |F|) array: for each permutation (row of S) the sorted list of image codes of F."""
        pw = np.int32(1) << S.astype(np.int32)                                    # (m, n)
        M = np.stack([pw[:, self.pts[m]].sum(axis=1) for m in F], axis=1)
        M.sort(axis=1)
        return M

    def is_canonical(self, F):
        Fs = np.array(sorted(F), dtype=np.int32)
        M = self.images(F, self.relevant_perms(F))
        return not rows_less(M, Fs).any()

    def automorphisms(self, F):
        """Aut(F) as an (m, n) int8 array of permutations (F must be canonical, or empty)."""
        if not F:
            return np.array(list(itertools.permutations(range(self.n))), dtype=np.int8)
        Fs = np.array(sorted(F), dtype=np.int32)
        S = self.relevant_perms(F)
        M = self.images(F, S)
        return S[np.all(M == Fs, axis=1)]

    def run(self):
        P = self.P

        def rec(F, last, cov, D):
            self.nodes += 1
            self.families.append((list(F), list(cov), D))
            for b in self.big:
                if b <= last:
                    continue
                if any(popcount(a & b) > 2 for a in F):
                    continue
                inc = comb(self.size[b] - 1, 2)
                if any(cov[p] + inc > P.cap_derived for p in self.pts[b]):
                    continue
                F2 = F + [b]
                if not self.is_canonical(F2):
                    continue
                cov2 = cov[:]
                for p in self.pts[b]:
                    cov2[p] += inc
                rec(F2, b, cov2, D + comb(self.size[b], 3) - 1)

        rec([], -1, [0] * self.n, 0)
        return self.families


# ----------------------------------------------------------------------------------------------------
# Lines: maximum number of pairwise <=1-intersecting blocks (3-blocks or rich blocks) under the caps
# ----------------------------------------------------------------------------------------------------
def line_sets(masks, P, need, collect=True, max_sets=100000):
    """masks: the rich blocks.  Returns (ell_max, list of line sets (as lists of masks) of size >= need)."""
    n = P.n
    covered = set()
    for m in masks:
        for t in itertools.combinations([i for i in range(n) if m >> i & 1], 3):
            covered.add(sum(1 << i for i in t))
    cands = []
    for t in itertools.combinations(range(n), 3):
        tm = sum(1 << i for i in t)
        if tm not in covered:
            cands.append((3, tm))
    for m in masks:
        cands.append((comb(popcount(m), 2), m))
    cands.sort(key=lambda x: (x[0], x[1]))
    N = len(cands)
    adj = [0] * N
    for i in range(N):
        for j in range(N):
            if i != j and popcount(cands[i][1] & cands[j][1]) <= 1:
                adj[i] |= 1 << j
    costs = [c for c, _ in cands]
    is3 = [c == 3 for c, _ in cands]
    best = [0]
    sets = []

    def ub(avail, cost, n3):
        cnt = popcount(avail)
        b = min(cnt, (P.cap_lines - cost) // 3 if P.cap_lines < INF else cnt)
        if P.max3lines < INF:
            rich = popcount(avail & ~tri_mask)
            b = min(b, (P.max3lines - n3) + rich)
        return b

    tri_mask = sum(1 << i for i in range(N) if is3[i])

    def rec(avail, chosen, cost, n3, threshold):
        k = len(chosen)
        if k > best[0]:
            best[0] = k
        if collect and k >= need:
            if len(sets) < max_sets:
                sets.append([cands[i][1] for i in chosen])
        a = avail
        while a:
            i = (a & -a).bit_length() - 1
            a &= a - 1
            ci = costs[i]
            if cost + ci > P.cap_lines or (is3[i] and n3 + 1 > P.max3lines):
                continue
            rest = avail & adj[i] & ~((1 << (i + 1)) - 1)
            thr = threshold() if callable(threshold) else threshold
            if k + 1 + ub(rest, cost + ci, n3 + is3[i]) < thr:
                continue
            rec(rest, chosen + [i], cost + ci, n3 + is3[i], threshold)

    full = (1 << N) - 1
    # pass 1: maximum (prune with best + 1)
    rec(full, [], 0, 0, lambda: best[0] + 1)
    ell = best[0]
    sets = []
    if collect and need <= ell:
        rec(full, [], 0, 0, max(need, 1))
    return ell, sets


# ----------------------------------------------------------------------------------------------------
# Phase 2: orderly extension of one canonical F_big by 4-blocks, canonicity under Aut(F_big)
# ----------------------------------------------------------------------------------------------------
class Phase2:
    def __init__(self, P, Fbig, cov0, D0, aut, verbose=False, max_leaves=INF):
        n = self.n = P.n
        self.P = P
        self.Fbig, self.D0 = list(Fbig), D0
        self.verbose = verbose
        self.subs = list(itertools.combinations(range(n), 4))                      # lex order = code order
        self.nc = len(self.subs)
        self.smask = np.array([sum(1 << i for i in s) for s in self.subs], dtype=np.int64)
        rank = np.full(1 << n, -1, dtype=np.int32)
        rank[self.smask] = np.arange(self.nc, dtype=np.int32)
        self.pts4 = np.array(self.subs, dtype=np.int64)                            # (nc, 4)
        self.rank = rank
        self.aut = aut
        self.m = len(aut)
        self.T = None                                                              # built lazily in run()
        smask = [int(x) for x in self.smask]
        self.allowed0 = np.array([all(popcount(smask[c] & b) <= 2 for b in Fbig) for c in range(self.nc)])
        self.compat = np.array([[popcount(a & b) <= 2 for b in smask] for a in smask])
        self.contains = np.array([[bool(m >> p & 1) for p in range(n)] for m in smask])   # (nc, n)
        self.triples = list(itertools.combinations(range(n), 3))
        trank = {t: i for i, t in enumerate(self.triples)}
        self.cand_triples = np.zeros((self.nc, len(self.triples)), dtype=bool)
        for c, s in enumerate(self.subs):
            for t in itertools.combinations(s, 3):
                self.cand_triples[c, trank[t]] = True
        self.pairs = list(itertools.combinations(range(n), 2))
        prank = {pr: i for i, pr in enumerate(self.pairs)}
        self.pt_inc = np.zeros((len(self.pairs), len(self.triples)), dtype=np.int64)
        for j, t in enumerate(self.triples):
            for pr in itertools.combinations(t, 2):
                self.pt_inc[prank[pr], j] = 1
        self.pair_pts = np.array(self.pairs, dtype=np.int64)
        self.cov = np.array(cov0, dtype=np.int64)
        self.deg4 = np.zeros(n, dtype=np.int64)
        for b in Fbig:                                   # a marked 4-block (--marked4 mode) counts for the t3 cap
            if popcount(b) == 4:
                for p in range(n):
                    if b >> p & 1:
                        self.deg4[p] += 1
        self.nodes = 0
        self.leaves = 0
        self.results = []
        self.max_leaves = max_leaves

    def root_bound(self):
        """Upper bound on the number of 4-blocks addable to F_big (used to skip hopeless families)."""
        P = self.P
        rem = np.minimum((P.cap_derived - self.cov) // 3, P.max4 - self.deg4)
        R = np.flatnonzero(self.allowed0)
        avail = self.contains[R].sum(axis=0)
        return int(np.minimum(rem, avail).sum() // 4)

    def build_T(self, chunk=20000):
        """T[sigma, c] = code of sigma(block c), built in chunks of permutations (memory-safe)."""
        T = np.empty((self.m, self.nc), dtype=np.int16)
        for s in range(0, self.m, chunk):
            a = self.aut[s:s + chunk].astype(np.int64)
            img = a[:, self.pts4]                                                  # (chunk, nc, 4)
            masks = (np.int64(1) << img).sum(axis=2)
            T[s:s + chunk] = self.rank[masks]
        assert (T >= 0).all()
        self.T = T

    def run(self):
        if self.m > 1 and self.allowed0.any():
            self.build_T()
        rel = np.zeros(self.m, dtype=bool) if self.m > 1 else None
        self.rec([], -1, self.D0, self.allowed0.copy(), rel)
        return self.results

    def rec(self, F4, last, D, allowed, rel):
        P = self.P
        self.nodes += 1
        if D >= P.D_min:
            self.leaf(F4, D)
        need = max(0, ceil_div(P.D_min - D, 3)) if P.D_min > -INF else 0
        rem = np.minimum((P.cap_derived - self.cov) // 3, P.max4 - self.deg4) if P.cap_derived < INF \
            else np.full(self.n, INF, dtype=np.int64)
        R = np.flatnonzero(allowed)
        R = R[R > last]
        if len(R) == 0:
            return
        sat = rem <= 0
        if sat.any():
            R = R[~self.contains[R][:, sat].any(axis=1)]
            if len(R) == 0:
                return
        if need > 0:
            # pair bound on rem
            tri_av = self.cand_triples[R].any(axis=0).astype(np.int64)
            thirds = self.pt_inc @ tri_av
            pair_cap = thirds // 2
            pcap = np.zeros(self.n, dtype=np.int64)
            np.add.at(pcap, self.pair_pts[:, 0], pair_cap)
            np.add.at(pcap, self.pair_pts[:, 1], pair_cap)
            rem = np.minimum(rem, pcap // 3)
            onehot = self.contains[R].astype(np.int64)
            suffix = np.cumsum(onehot[::-1], axis=0)[::-1]
            bound = np.minimum(suffix, rem).sum(axis=1) // 4
            if bound[0] < need:
                return
        else:
            bound = None
        c0 = F4[0] if F4 else None
        for j in range(len(R)):
            if bound is not None and bound[j] < need:
                break
            c = int(R[j])
            newF = F4 + [c]
            if self.m > 1:
                cc0 = c0 if c0 is not None else c
                newrel = rel | (self.T[:, c] <= cc0)
                idx = np.flatnonzero(newrel)
                rows = self.T[idx[:, None], np.array(newF, dtype=np.int64)[None, :]]
                rows.sort(axis=1)
                if rows_less(rows, np.array(newF, dtype=np.int16)).any():
                    continue
            else:
                newrel = None
            pts = list(self.subs[c])
            self.cov[pts] += 3
            self.deg4[pts] += 1
            self.rec(newF, c, D + 3, allowed & self.compat[c], newrel)
            self.cov[pts] -= 3
            self.deg4[pts] -= 1

    def leaf(self, F4, D):
        P = self.P
        self.leaves += 1
        if self.leaves > self.max_leaves:
            return
        if P.ell_max >= INF:
            self.results.append({"Fbig": list(self.Fbig), "F4": [int(self.smask[c]) for c in F4], "D": D})
            return
        need_l = P.N3 - P.target - D
        masks = list(self.Fbig) + [int(self.smask[c]) for c in F4]
        ell, sets = line_sets(masks, P, need_l)
        if ell >= need_l:
            n = self.n
            dec = lambda m: [i for i in range(n) if m >> i & 1]
            rec = {"n": n, "mode": P.mode, "blocks": [dec(m) for m in masks],
                   "sizes": sorted([popcount(m) for m in masks], reverse=True), "D": D, "ell_max": ell,
                   "count_min": P.N3 - D - ell,
                   "degrees": [sum(1 for m in masks if m >> p & 1) for p in range(n)],
                   "line_sets": [[dec(l) for l in L] for L in sets]}
            self.results.append(rec)
            if self.verbose:
                print(f"    CANDIDATE sizes={rec['sizes']} D={D} ell_max={ell} count={rec['count_min']} "
                      f"degrees={rec['degrees']} #line-sets>= {need_l}: {len(sets)}", flush=True)


# ----------------------------------------------------------------------------------------------------
def run(P, include_empty=False, family_filter=None, verbose=True, log=print, out=None):
    ph1 = Phase1(P)
    fams = ph1.run()
    log(f"phase 1: {len(fams)} canonical big-block families (incl. empty)")
    results = []
    stats = []
    total_nodes = 0
    for idx, (Fbig, cov, D0) in enumerate(fams):
        if not Fbig and not include_empty:
            continue
        if family_filter is not None and idx not in family_filter:
            continue
        aut = ph1.automorphisms(Fbig)
        ph2 = Phase2(P, Fbig, cov, D0, aut, verbose=verbose)
        rb = ph2.root_bound()
        sizes = [ph1.size[m] for m in Fbig]
        if P.D_min > -INF and D0 + 3 * rb < P.D_min:
            log(f"[{idx}] F_big sizes={sizes} blocks={[ph1.pts[m] for m in Fbig]} D0={D0} |Aut|={len(aut)} "
                f"root bound {rb} four-blocks -> D <= {D0 + 3 * rb} < D_min={P.D_min}: skipped")
            stats.append({"idx": idx, "sizes": sizes, "D0": D0, "aut": int(len(aut)), "root_bound": rb, "nodes": 0,
                          "leaves": 0, "skipped": True})
            continue
        t0 = time.time()
        log(f"[{idx}] F_big sizes={sizes} blocks={[ph1.pts[m] for m in Fbig]} D0={D0} |Aut|={len(aut)} "
            f"root bound {rb} (need >= {max(0, ceil_div(P.D_min - D0, 3)) if P.D_min > -INF else 0}) ...")
        res = ph2.run()
        total_nodes += ph2.nodes
        log(f"      nodes={ph2.nodes} leaves={ph2.leaves} candidates={len(res)} time={time.time() - t0:.1f}s")
        stats.append({"idx": idx, "sizes": sizes, "D0": D0, "aut": int(len(aut)), "root_bound": rb,
                      "nodes": ph2.nodes, "leaves": ph2.leaves, "candidates": len(res),
                      "time": round(time.time() - t0, 1), "skipped": False})
        results.extend(res)
        if out:
            json.dump({"params": P.describe(), "results": results, "stats": stats, "complete": False}, open(out, "w"), indent=1)
    log(f"total phase-2 nodes {total_nodes}; candidate structures with count <= {P.target}: {len(results)}")
    if out:
        json.dump({"params": P.describe(), "results": results, "stats": stats, "complete": True}, open(out, "w"), indent=1)
    return results, stats, fams


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("target", type=int)
    ap.add_argument("--mode", default="table")
    ap.add_argument("--big", default=None)
    ap.add_argument("--include-empty", action="store_true")
    ap.add_argument("--families", default=None, help="comma-separated phase-1 family indices to run")
    ap.add_argument("--out", default=None)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--list", action="store_true", help="phase 1 + root bounds only")
    ap.add_argument("--marked4", action="store_true", help="all-4-block case with a marked first block")
    a = ap.parse_args()
    big = [int(x) for x in a.big.split(",")] if a.big else None
    P = Params(a.n, a.target, a.mode, big)
    print(P.describe(), flush=True)
    ff = set(int(x) for x in a.families.split(",")) if a.families else None
    if a.marked4:
        # all-4-block case: mark the lexicographically first block {0,1,2,3}; canonicity of the rest under its
        # stabiliser S_4 x S_6.  Every all-4-block family is visited at least once (once per Aut-orbit of blocks).
        ph1 = Phase1(P)
        b0 = (1 << 4) - 1
        aut = ph1.block_perms(b0)
        cov0 = [3 if p < 4 else 0 for p in range(P.n)]
        ph2 = Phase2(P, [b0], cov0, 3, aut, verbose=not a.quiet)
        rb = ph2.root_bound()
        print(f"marked first 4-block {ph1.pts[b0]}: |Stab|={len(aut)} root bound {rb} need >= "
              f"{max(0, ceil_div(P.D_min - 3, 3)) if P.D_min > -INF else 0} #cands={int(ph2.allowed0.sum())}", flush=True)
        t0 = time.time()
        res = ph2.run()
        print(f"nodes={ph2.nodes} leaves={ph2.leaves} candidates={len(res)} time={time.time() - t0:.1f}s", flush=True)
        if a.out:
            json.dump({"params": P.describe() + " marked4", "results": res, "complete": True}, open(a.out, "w"), indent=1)
        return
    if a.list:
        ph1 = Phase1(P)
        fams = ph1.run()
        print(f"phase 1: {len(fams)} canonical big-block families (incl. empty)")
        for idx, (Fbig, cov, D0) in enumerate(fams):
            if not Fbig:
                print(f"[{idx}] empty family (Aut = S_n)"); continue
            aut = ph1.automorphisms(Fbig)
            ph2 = Phase2(P, Fbig, cov, D0, aut)
            rb = ph2.root_bound()
            need = max(0, ceil_div(P.D_min - D0, 3)) if P.D_min > -INF else 0
            print(f"[{idx}] sizes={[ph1.size[m] for m in Fbig]} blocks={[ph1.pts[m] for m in Fbig]} D0={D0} "
                  f"|Aut|={len(aut)} root_bound={rb} need={need} {'SKIP' if rb < need else 'RUN'} "
                  f"#cands={int(ph2.allowed0.sum())}")
        return
    results, stats, fams = run(P, include_empty=a.include_empty, family_filter=ff, verbose=not a.quiet,
                               log=lambda s: print(s, flush=True), out=a.out)


if __name__ == "__main__":
    main()
