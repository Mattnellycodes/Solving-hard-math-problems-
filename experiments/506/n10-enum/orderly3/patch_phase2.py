#!/usr/bin/env python3
"""One-off patch: replace class Phase2 in oe.py by the version with suffix/thirds bounds,
child-bound-before-canonicity and incremental relevant automorphisms."""
src = open('oe.py').read()
start = src.index("class Phase2:")
end = src.index("# ----------------------------------------------------------------------------------------------\ndef run_enum")
new_class = '''class Phase2:
    """Orderly extension of a canonical big-block family F_big by 4-blocks.

    Codes of 4-blocks = their index in the list of 4-subset bitmasks sorted increasingly.
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
        self.masks4 = [list_to_mask(c) for c in itertools.combinations(range(n), 4)]
        self.masks4.sort()
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


'''
src = src[:start] + new_class + src[end:]
open('oe.py', 'w').write(src)
print("Phase2 rewritten")
