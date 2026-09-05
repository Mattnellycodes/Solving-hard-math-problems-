#!/usr/bin/env python3
"""
all4.py -- orderly enumeration of families consisting of 4-blocks only (largest block 4), with a
full-S_n canonicity test that is cheap enough for n = 10.

Codes: 4-blocks are numbered in lexicographic order of their sorted point tuples (code 0 =
{0,1,2,3}; all blocks through point 0 come first).  Canonical form of a family = the
lexicographically smallest sorted tuple of codes over S_n.  In a canonical family the first code
is 0.  Test for a family F' = (0 = F'_0, F'_1, ...):
a permutation sigma with sorted(sigma(F')) <= F' must map some block c to {0,1,2,3} and some other
block c' to a block m with code(m) <= F'_1 and |m & {0,1,2,3}| = |c & c'| (the second smallest
image).  For each such (c, c', m) the permutations are tau o rho with rho a fixed bijection
(c -> {0,1,2,3}, c' -> m, respecting the intersection) and tau in the ordered stabiliser of
({0,1,2,3}, m).  All of them are enumerated, so the test is exactly the brute-force test (validated
against brute force for n = 7, 8 in validate_all4()).  Removing the largest code of a canonical
family leaves a canonical family, so the DFS adding blocks in increasing code order and keeping
canonical children visits each isomorphism class exactly once.

Bounds (as in oe.Phase2): per point a_p = min((capD-used_p)//3, max4-d4_p, #cands through p,
floor(sum_q floor(thirds_pq/2)/3)); #additional <= sum_p a_p // 4; suffix bound over the candidate
list; the child's bound is checked before its (expensive) canonicity test.
"""
import sys, json, time, itertools
from math import comb
import numpy as np
from oe import (Caps, popcount, mask_to_list, list_to_mask, stab_perms, images, all_perms,
                canon_bruteforce, line_sets)


class All4:
    def __init__(self, n, caps, need_total, log=None):
        self.n, self.caps, self.need_total = n, caps, need_total
        self.masks = [list_to_mask(c) for c in itertools.combinations(range(n), 4)]   # lex order
        self.M = len(self.masks)
        self.mask_arr = np.array(self.masks, dtype=np.int32)
        self.pts = np.array([mask_to_list(m, n) for m in self.masks], dtype=np.int8)
        self.onehot = np.zeros((self.M, n), dtype=np.int16)
        for j in range(self.M):
            self.onehot[j, self.pts[j]] = 1
        self.code_of = np.full(1 << n, -1, dtype=np.int32)
        self.code_of[self.mask_arr] = np.arange(self.M, dtype=np.int32)
        self.pc = np.array([popcount(m) for m in range(1 << n)], dtype=np.int8)
        # triples / pairs for the thirds bound
        self.triples = [list_to_mask(c) for c in itertools.combinations(range(n), 3)]
        tri_index = {m: i for i, m in enumerate(self.triples)}
        self.tri_of_block = np.array([[tri_index[m & ~(1 << p)] for p in mask_to_list(m, n)]
                                      for m in self.masks], dtype=np.int32)
        pairs = list(itertools.combinations(range(n), 2))
        self.pair_tri = np.array([[tri_index[(1 << a) | (1 << b) | (1 << x)] for x in range(n) if x != a and x != b]
                                  for (a, b) in pairs], dtype=np.int32)
        self.pt_pair = np.zeros((n, len(pairs)), dtype=np.int64)
        for k, (a, b) in enumerate(pairs):
            self.pt_pair[a, k] = 1
            self.pt_pair[b, k] = 1
        # stabiliser of {0,1,2,3} and its code-image table
        self.stab4 = stab_perms(n, 4)
        self.stab_img = self.code_of[images(self.stab4, self.masks)]        # (S, M)
        # candidate second blocks by intersection size with {0,1,2,3}, in code order
        self.Mi = {i: [(self.code_of[m], m) for m in self.masks if m != 15 and popcount(m & 15) == i]
                   for i in (0, 1, 2)}
        self.sub_cache = {}
        self.Dmin = need_total - caps.ellmax
        self.nodes = self.evals = self.tests = 0
        self.results = []
        self.log = log

    def sub_stab(self, m):
        if m not in self.sub_cache:
            j = self.code_of[m]
            sel = self.stab_img[:, j] == j
            self.sub_cache[m] = self.stab_img[sel]
        return self.sub_cache[m]

    def rho(self, c, cp, m):
        """permutation sending block c -> {0,1,2,3} and c' -> m (masks), as a list."""
        n = self.n
        inter = [p for p in range(n) if (c >> p & 1) and (cp >> p & 1)]
        c_only = [p for p in range(n) if (c >> p & 1) and not (cp >> p & 1)]
        cp_only = [p for p in range(n) if (cp >> p & 1) and not (c >> p & 1)]
        rest = [p for p in range(n) if not (c >> p & 1) and not (cp >> p & 1)]
        t_inter = [p for p in range(4) if m >> p & 1]
        t_c = [p for p in range(4) if not (m >> p & 1)]
        t_cp = [p for p in range(4, n) if m >> p & 1]
        t_rest = [p for p in range(4, n) if not (m >> p & 1)]
        r = [0] * n
        for src, dst in ((inter, t_inter), (c_only, t_c), (cp_only, t_cp), (rest, t_rest)):
            assert len(src) == len(dst)
            for a, b in zip(src, dst):
                r[a] = b
        return r

    def canonical(self, F):
        """F: sorted list of codes."""
        k = len(F)
        if F[0] != 0:
            return False
        if k == 1:
            return True
        self.tests += 1
        F1 = F[1]
        target = np.array(F, dtype=np.int32)
        Fm = [self.masks[c] for c in F]
        for i in (2, 1, 0):
            ms = [m for (code, m) in self.Mi[i] if code <= F1]
            if not ms:
                continue
            pairs = [(c, cp) for c in Fm for cp in Fm if c != cp and popcount(c & cp) == i]
            if not pairs:
                continue
            for m in ms:
                T = self.sub_stab(m)                                          # (S_m, M)
                rho_arr = np.array([self.rho(c, cp, m) for (c, cp) in pairs], dtype=np.int8)
                fam_idx = np.sort(self.code_of[images(rho_arr, Fm)], axis=1)   # (npairs, k)
                rows = np.sort(T[:, fam_idx], axis=2)                         # (S_m, npairs, k)
                lt = np.zeros(rows.shape[:2], dtype=bool)
                eq = np.ones(rows.shape[:2], dtype=bool)
                for j in range(k):
                    col = rows[:, :, j]
                    lt |= eq & (col < target[j])
                    eq &= (col == target[j])
                if lt.any():
                    return False
        return True

    # ----- bounds
    def point_caps(self, cands, used, d4):
        caps = self.caps
        a = np.minimum((caps.capD - used) // 3, caps.max4 - d4).astype(np.int64)
        if len(cands) == 0:
            return np.zeros(self.n, dtype=np.int64)
        tri_av = np.zeros(len(self.triples), dtype=np.int64)
        tri_av[self.tri_of_block[cands].ravel()] = 1
        thirds = tri_av[self.pair_tri].sum(axis=1)
        per_pt = self.pt_pair @ (thirds // 2)
        a = np.minimum(a, per_pt // 3)
        a = np.minimum(a, self.onehot[cands].sum(axis=0))
        return np.maximum(a, 0)

    def bound_extra(self, cands, used, d4):
        if len(cands) == 0:
            return 0
        return min(len(cands), int(self.point_caps(cands, used, d4).sum()) // 4)

    def evaluate(self, F):
        self.evals += 1
        Fm = [self.masks[c] for c in F]
        D = 3 * len(F)
        need = self.need_total - D
        lmax, sets = line_sets(self.n, Fm, self.caps, need)
        if lmax >= need:
            cnt = comb(self.n, 3) - D - lmax
            deg = [sum(1 for B in Fm if B >> p & 1) for p in range(self.n)]
            rec = {"n": self.n, "mode": self.caps.mode, "blocks": [mask_to_list(B, self.n) for B in Fm],
                   "sizes": [4] * len(F), "D": D, "ell_max": lmax, "count_min": cnt, "degrees": deg,
                   "line_sets": [[mask_to_list(L, self.n) for L in S] for S in sets]}
            self.results.append(rec)
            if self.log:
                print(f"    CANDIDATE: b4={len(F)} D={D} ell_max={lmax} count={cnt} degrees={deg} "
                      f"#line_sets={len(sets)}", file=self.log, flush=True)

    def run(self):
        n = self.n
        used = np.zeros(n, dtype=np.int32)
        d4 = np.zeros(n, dtype=np.int32)
        self.t0 = time.time()
        self._rec([], np.arange(self.M), used, d4)
        return self.results

    def _rec(self, F, cands, used, d4):
        caps = self.caps
        self.nodes += 1
        if self.log and self.nodes % 20000 == 0:
            print(f"    ... nodes={self.nodes} tests={self.tests} evals={self.evals} depth={len(F)} "
                  f"[{time.time()-self.t0:.0f}s]", file=self.log, flush=True)
        D = 3 * len(F)
        if D >= self.Dmin:
            self.evaluate(F)
        need_extra = max(0, -(-(self.need_total - caps.ellmax - D) // 3))
        if len(cands) == 0:
            return
        pts, pc, mask_arr = self.pts, self.pc, self.mask_arr
        if need_extra > 0:
            a = self.point_caps(cands, used, d4)
            suffix = np.cumsum(self.onehot[cands][::-1], axis=0)[::-1]
            sbound = np.minimum(suffix, a).sum(axis=1) // 4
        for pos in range(len(cands)):
            if need_extra > 0 and sbound[pos] < need_extra:
                break
            c = int(cands[pos])
            used2 = used.copy()
            d42 = d4.copy()
            for p in pts[c]:
                used2[p] += 3
                d42[p] += 1
            point_ok = (used2 + 3 <= caps.capD) & (d42 < caps.max4)
            rest = cands[pos + 1:]
            keep = (pc[mask_arr[rest] & mask_arr[c]] <= 2) & point_ok[pts[rest]].all(axis=1)
            rest = rest[keep]
            if need_extra > 1 and self.bound_extra(rest, used2, d42) < need_extra - 1:
                continue
            child = F + [c]
            if not self.canonical(child):
                continue
            self._rec(child, rest, used2, d42)


def canon_bruteforce_codes(a4, fam_masks, perms):
    """lexicographically smallest sorted code tuple over all of S_n (validation only)."""
    rows = np.sort(a4.code_of[images(perms, list(fam_masks))], axis=1)
    order = np.lexsort(rows.T[::-1])
    return [int(x) for x in rows[order[0]]]


def validate_all4():
    """canonical() must agree with brute force (in code order) on random 4-block families for
    n = 7, 8; and the DFS must produce exactly the brute-force isomorphism classes (n = 7)."""
    import random
    random.seed(7)
    for n in (7, 8):
        perms = all_perms(n)
        a4 = All4(n, Caps(n, 'none'), 0)
        masks = list(a4.masks)
        nt = 0
        for trial in range(400 if n == 7 else 200):
            random.shuffle(masks)
            fam = []
            for b in masks[:random.randint(1, 12)]:
                if all(popcount(b & B) <= 2 for B in fam):
                    fam.append(b)
            codes = sorted(int(a4.code_of[b]) for b in fam)
            cf = canon_bruteforce_codes(a4, fam, perms)
            for f in (codes, cf):
                truth = (f == cf)
                got = a4.canonical(f)
                if truth != got:
                    print("MISMATCH", n, f, cf, truth, got)
                    return False
                nt += 1
        print(f"validate_all4 n={n}: {nt} canonicity tests agree with brute force", flush=True)
    n = 7
    perms = all_perms(n)
    masks = sorted(list_to_mask(c) for c in itertools.combinations(range(n), 4))
    fams = []

    def rec(start, fam):
        fams.append(tuple(fam))
        for i in range(start, len(masks)):
            if all(popcount(masks[i] & B) <= 2 for B in fam):
                fam.append(masks[i]); rec(i + 1, fam); fam.pop()
    rec(0, [])
    classes = set(canon_bruteforce(n, f, perms) for f in fams)
    a4 = All4(n, Caps(n, 'none'), 0)
    visited = []
    a4.evaluate = lambda F: visited.append(tuple(sorted(a4.masks[c] for c in F)))
    a4.Dmin = -10 ** 9
    a4.caps.ellmax = 10 ** 9
    a4.run()
    vs = set(canon_bruteforce(n, f, perms) for f in visited)
    ok = (vs == classes) and len(visited) == len(vs)
    print(f"validate_all4 n=7: labelled 4-block families={len(fams)}, brute-force classes={len(classes)}, "
          f"orderly visited={len(visited)} distinct classes={len(vs)} => {'MATCH' if ok else 'MISMATCH'}", flush=True)
    return ok


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('n', type=int)
    ap.add_argument('target', type=int)
    ap.add_argument('--mode', default='sg')
    ap.add_argument('--out', default=None)
    ap.add_argument('--validate', action='store_true')
    a = ap.parse_args()
    if a.validate:
        print("VALIDATION", "OK" if validate_all4() else "FAILED")
        sys.exit(0)
    caps = Caps(a.n, a.mode)
    need_total = comb(a.n, 3) - a.target
    print(f"all-4-block case n={a.n} target={a.target} mode={a.mode} need D+ell>={need_total} {caps}", flush=True)
    t0 = time.time()
    a4 = All4(a.n, caps, need_total, log=sys.stdout)
    res = a4.run()
    print(f"nodes={a4.nodes} canonicity tests={a4.tests} evals={a4.evals} candidates={len(res)} "
          f"[{time.time()-t0:.1f}s]", flush=True)
    for i, r in enumerate(res):
        print(f"  cand {i}: b4={len(r['blocks'])} D={r['D']} ell_max={r['ell_max']} count={r['count_min']} "
              f"degrees={r['degrees']}\n     blocks={r['blocks']}")
    if a.out:
        json.dump({"n": a.n, "target": a.target, "mode": a.mode, "caps": str(caps), "sizes": [],
                   "results": res}, open(a.out, 'w'), indent=1)
