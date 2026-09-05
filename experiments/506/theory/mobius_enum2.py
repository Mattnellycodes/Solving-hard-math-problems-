"""Two-phase orderly enumeration of Mobius block structures (for n = 9, 10 where the one-phase
enumerator of mobius_enum.py is too slow).

Phase 1: canonical (lex-min under S_n) families F_big of blocks with sizes in BIG (default 5..n-1),
pairwise sharing <= 2 points, subject to the derived Sylvester-Gallai cap (C2).
Phase 2: for each F_big, orderly generation of 4-block extensions F4.  The pair (F_big, F4) is
canonical iff the concatenated sequence (sorted F_big, sorted F4) is lex-minimal under S_n; because
F_big is already lex-minimal this only involves permutations in Aut(F_big), which we compute once.
The parent of a canonical pair (remove the largest 4-block) is canonical (proof in REPORT.md, App. A),
so orderly generation is exhaustive: every isomorphism class of (F_big, F4) is visited exactly once.

The case F_big = empty (all blocks of size 4) is handled separately (flag --empty), because
Aut(empty) = S_n makes the canonicity test slow; for n = 9 it is excluded by theory (REPORT Prop 5.2,
which uses the non-realisability of the (8_3) Mobius-Kantor configuration).

Usage: python3 mobius_enum2.py n target [--big 5,6,7,8] [--empty] [--out file.json]
"""
import sys, json, itertools, argparse, time
from math import comb
import numpy as np
from mobius_enum import o_lower, popcount, Enumerator


class Enumerator2:
    def __init__(self, n, target, big_sizes, verbose=True, o_func=o_lower):
        self.n = n; self.target = target; self.verbose = verbose
        self.N3 = comb(n, 3)
        self.cap_derived = comb(n - 1, 2) - o_func(n - 1)
        self.cap_lines = comb(n, 2) - o_func(n)
        self.ell_max = self.cap_lines // 3
        self.D_min = self.N3 - target - self.ell_max
        self.big = [c for c in range(1 << n) if popcount(c) in big_sizes]
        self.four = [c for c in range(1 << n) if popcount(c) == 4]
        self.size = {c: popcount(c) for c in range(1 << n)}
        self.deficit = {c: comb(self.size[c], 3) - 1 for c in range(1 << n)}
        self.pair_index = {}
        # permutations as arrays; Pow[i, j] = 2^{perm_j(i)}
        perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64)   # (P, n)
        self.P = len(perms)
        self.Pow = (1 << perms).T.astype(np.int64)                                # (n, P)
        codes = np.arange(1 << n, dtype=np.int64)
        self.Bits = ((codes[:, None] >> np.arange(n)) & 1).astype(np.int64)     # (2^n, n)
        self.results = []
        self.nodes = 0
        self.helper = Enumerator.__new__(Enumerator)   # reuse line_sets / uncovered_triples
        self.helper.n = n; self.helper.cap_lines = self.cap_lines
        self.helper.decode = lambda code: [i for i in range(n) if code >> i & 1]

    def images(self, F, cols=None):
        """(k, P) array of image codes of the blocks of F under all permutations (or given columns)."""
        B = self.Bits[np.array(F, dtype=np.int64)]                                 # (k, n)
        Pw = self.Pow if cols is None else self.Pow[:, cols]
        return B @ Pw

    def is_canonical_full(self, F):
        Fs = sorted(F)
        M = np.sort(self.images(Fs), axis=0)                                       # (k, P)
        alive = np.ones(M.shape[1], dtype=bool)
        for j, f in enumerate(Fs):
            if np.any(alive & (M[j] < f)):
                return False
            alive &= (M[j] == f)
        return True

    def automorphisms(self, F):
        """indices of permutations fixing the family F (as a set)."""
        if not F:
            return np.arange(self.P)
        Fs = sorted(F)
        M = np.sort(self.images(Fs), axis=0)
        ok = np.all(M == np.array(Fs)[:, None], axis=0)
        return np.where(ok)[0]

    def is_canonical_rel(self, F4, aut):
        """canonicity of (F_big, F4) given Aut(F_big) column indices."""
        if len(aut) == 1:
            return True
        Fs = sorted(F4)
        M = np.sort(self.images(Fs, aut), axis=0)
        alive = np.ones(M.shape[1], dtype=bool)
        for j, f in enumerate(Fs):
            if np.any(alive & (M[j] < f)):
                return False
            alive &= (M[j] == f)
        return True

    # ---------------- phase 1 ----------------
    def phase1(self, include_empty):
        out = []
        def rec(F, last, cov, D):
            if F or include_empty:
                out.append((list(F), list(cov), D))
            for b in self.big:
                if b <= last: continue
                if any(popcount(a & b) > 2 for a in F): continue
                k = self.size[b]; inc = comb(k - 1, 2)
                if any(cov[p] + inc > self.cap_derived for p in range(self.n) if b >> p & 1): continue
                F2 = F + [b]
                if not self.is_canonical_full(F2): continue
                cov2 = cov[:]
                for p in range(self.n):
                    if b >> p & 1: cov2[p] += inc
                rec(F2, b, cov2, D + self.deficit[b])
        rec([], -1, [0] * self.n, 0)
        return out

    # ---------------- phase 2 ----------------
    def phase2(self, Fbig, cov0, D0):
        aut = self.automorphisms(Fbig)
        if self.verbose:
            print(f"F_big sizes={[self.size[b] for b in Fbig]} D0={D0} |Aut|={len(aut)}", flush=True)
        cands4 = [b for b in self.four if all(popcount(a & b) <= 2 for a in Fbig)]
        cov = cov0[:]

        def rec(F4, last, D):
            self.nodes += 1
            if D >= self.D_min:
                self.process(Fbig + F4, D)
            # candidates
            ext = []
            for b in cands4:
                if b <= last: continue
                if any(popcount(a & b) > 2 for a in F4): continue
                if any(cov[p] + 3 > self.cap_derived for p in range(self.n) if b >> p & 1): continue
                ext.append(b)
            # potential: number of further 4-blocks
            cap_blocks = sum((self.cap_derived - cov[p]) // 3 for p in range(self.n)) // 4
            pot = min(len(ext), cap_blocks)
            if D + 3 * pot < self.D_min:
                return
            for idx, b in enumerate(ext):
                if D + 3 * min(len(ext) - idx, cap_blocks) < self.D_min:
                    break
                F2 = F4 + [b]
                if not self.is_canonical_rel(F2, aut):
                    continue
                for p in range(self.n):
                    if b >> p & 1: cov[p] += 3
                rec(F2, b, D + 3)
                for p in range(self.n):
                    if b >> p & 1: cov[p] -= 3
        rec([], -1, D0)

    def process(self, F, D):
        need = max(self.N3 - self.target - D, 0)
        ell, sets = self.helper.line_sets(F, need)
        if ell < need:
            return
        count = self.N3 - D - ell
        rec = {"n": self.n, "blocks": [sorted(self.helper.decode(b)) for b in sorted(F)],
               "sizes": sorted([self.size[b] for b in F], reverse=True), "D": D,
               "ell_max": ell, "count_min": count,
               "degrees": [sum(1 for b in F if b >> p & 1) for p in range(self.n)],
               "line_sets": [[sorted(self.helper.decode(l)) for l in L] for L in sets]}
        self.results.append(rec)
        if self.verbose:
            print(f"  candidate: sizes={rec['sizes']} D={D} ell_max={ell} count={count} "
                  f"degrees={rec['degrees']} #line-sets={len(sets)}", flush=True)

    def run(self, include_empty=False):
        fams = self.phase1(include_empty)
        if self.verbose:
            print(f"phase 1: {len(fams)} canonical big-block families", flush=True)
        for Fbig, cov, D0 in fams:
            # quick prune: max possible deficit with 4-blocks
            cap_blocks = sum((self.cap_derived - cov[p]) // 3 for p in range(self.n)) // 4
            if D0 + 3 * cap_blocks < self.D_min:
                continue
            t0 = time.time()
            self.phase2(Fbig, cov, D0)
            if self.verbose:
                print(f"   done in {time.time()-t0:.1f}s, nodes so far {self.nodes}", flush=True)
        return self.results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int); ap.add_argument("target", type=int)
    ap.add_argument("--big", default=None); ap.add_argument("--empty", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    big = [int(x) for x in a.big.split(",")] if a.big else list(range(5, a.n))
    E = Enumerator2(a.n, a.target, big)
    print(f"n={a.n} target={a.target}: C(n,3)={E.N3} cap_derived={E.cap_derived} cap_lines={E.cap_lines} "
          f"ell_max={E.ell_max} D_min={E.D_min} big sizes={big} include_empty={a.empty}", flush=True)
    res = E.run(include_empty=a.empty)
    print(f"nodes {E.nodes}; candidates with count <= {a.target}: {len(res)}")
    if a.out:
        json.dump(res, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
