"""Independent orderly enumeration (Read-Faradzev) of abstract block structures for Erdos #506.
Family F = blocks of size 4..n-1 (bitmasks), pairwise sharing <= 2 points, per-point pair coverage
sum_{B∋p} C(|B|-1,2) <= capd.  Canonical form = lexicographically least sorted code tuple under S_n.
Soundness of orderly generation: if F is canonical then F minus its largest code is canonical
(inserting one element into a sorted list never increases the value at any index); hence every
isomorphism class is reached exactly once by extending canonical families with larger codes.
For each family with D >= D_min the exact maximum number of lines ell_max(F) is computed by
branch and bound over L ⊆ F ∪ {uncovered triples}, pairwise sharing <= 1 point, sum C(|l|,2) <= capl.
Usage: python3 orderly_audit.py n target MODE [--validate]     MODE in strong|weak|none
"""
import sys, itertools, json, time
from math import comb
import numpy as np

# rigorous lower bounds for o(m) = min number of ordinary lines of m non-collinear points
def o_val(m, mode):
    if mode == "none": return 0
    if mode == "weak": return 1                      # Sylvester-Gallai
    table = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4}      # classical values (Kelly-Moser gives >= 3m/7)
    return table.get(m, max(1, -(-3 * m // 7)))

pc = lambda x: bin(x).count("1")

class Orderly:
    def __init__(self, n, target, mode):
        self.n, self.target = n, target
        self.capd = comb(n - 1, 2) - o_val(n - 1, mode)
        self.capl = comb(n, 2) - o_val(n, mode)
        self.ellmax_global = self.capl // 3
        self.Dmin = comb(n, 3) - target - self.ellmax_global
        self.cands = [c for c in range(1 << n) if 4 <= pc(c) <= n - 1]
        self.deficit = {c: comb(pc(c), 3) - 1 for c in self.cands}
        self.cov1 = {c: comb(pc(c) - 1, 2) for c in self.cands}    # pairs covered at each point of c
        perms = list(itertools.permutations(range(n)))
        codes = np.arange(1 << n)
        bits = (codes[:, None] >> np.arange(n)) & 1                 # (2^n, n)
        img = np.zeros((len(perms), n), dtype=np.int64)
        for j, s in enumerate(perms):
            img[j] = [1 << s[i] for i in range(n)]
        self.PI = (bits @ img.T).T.astype(np.int32)                # PI[j, code] = image of code under perm j
        self.nodes = 0; self.found = []; self.families_checked = 0

    def canonical(self, F):                       # F sorted list of codes
        Fs = np.array(F, dtype=np.int32)
        M = np.sort(self.PI[:, Fs], axis=1)
        eq = np.ones(M.shape[0], dtype=bool)
        for j in range(len(F)):
            col = M[:, j]
            if np.any(eq & (col < Fs[j])):
                return False
            eq &= (col == Fs[j])
        return True

    def run(self):
        self.dfs([], -1, [0] * self.n, 0)
        return self.found

    def dfs(self, F, last, cov, D):
        self.nodes += 1
        if D >= self.Dmin:
            self.check(F, D)
        budget = sum(self.capd - c for c in cov)
        if D + budget // 3 < self.Dmin:
            return
        ext = []
        for c in self.cands:
            if c <= last: continue
            if any(pc(c & b) > 2 for b in F): continue
            if any(cov[p] + self.cov1[c] > self.capd for p in range(self.n) if c >> p & 1): continue
            ext.append(c)
        pot = sum(self.deficit[c] for c in ext)
        for i, c in enumerate(ext):
            if D + pot < self.Dmin:
                break
            pot -= self.deficit[c]
            F2 = sorted(F + [c])
            if not self.canonical(F2):
                continue
            for p in range(self.n):
                if c >> p & 1: cov[p] += self.cov1[c]
            self.dfs(F + [c], c, cov, D + self.deficit[c])
            for p in range(self.n):
                if c >> p & 1: cov[p] -= self.cov1[c]

    def ell_max(self, F):
        n = self.n
        unc = [sum(1 << i for i in t) for t in itertools.combinations(range(n), 3)
               if not any((sum(1 << i for i in t) & b) == sum(1 << i for i in t) for b in F)]
        cand = sorted(F) + unc
        w = [comb(pc(c), 2) for c in cand]
        ok = [[pc(cand[i] & cand[j]) <= 1 for j in range(len(cand))] for i in range(len(cand))]
        best = [0]; sets = []
        def rec(i, chosen, used):
            if len(chosen) > best[0]:
                best[0] = len(chosen); sets.clear()
            if len(chosen) == best[0]:
                sets.append([cand[k] for k in chosen])
            avail = [j for j in range(i, len(cand)) if all(ok[k][j] for k in chosen) and used + w[j] <= self.capl]
            if len(chosen) + len(avail) < best[0]:
                return
            for idx, j in enumerate(avail):
                if len(chosen) + len(avail) - idx < best[0]:
                    break
                rec(j + 1, chosen + [j], used + w[j])
        rec(0, [], 0)
        return best[0], sets

    def check(self, F, D):
        self.families_checked += 1
        ell, sets = self.ell_max(F)
        count = comb(self.n, 3) - D - ell
        if count <= self.target:
            dec = lambda c: [i for i in range(self.n) if c >> i & 1]
            rec = {"blocks": [dec(b) for b in sorted(F)], "sizes": sorted((pc(b) for b in F), reverse=True),
                   "D": D, "ell_max": ell, "count": count,
                   "degrees": [sum(1 for b in F if b >> p & 1) for p in range(self.n)],
                   "num_max_line_sets": len(sets), "max_line_sets": [[dec(l) for l in L] for L in sets]}
            self.found.append(rec)
            print(f"  FOUND sizes={rec['sizes']} D={D} ell_max={ell} count={count} degrees={rec['degrees']} "
                  f"#max-line-sets={len(sets)}", flush=True)

def brute_classes(n, maxk):
    """all isomorphism classes of families of >=4-blocks (pairwise <=2), labelled enumeration + canonical form"""
    cands = [c for c in range(1 << n) if 4 <= pc(c) <= n - 1]
    perms = list(itertools.permutations(range(n)))
    def img(c, s): return sum(1 << s[i] for i in range(n) if c >> i & 1)
    def canon(F): return min(tuple(sorted(img(c, s) for c in F)) for s in perms)
    classes = set()
    def rec(F, last):
        classes.add(canon(F))
        if len(F) >= maxk: return
        for c in cands:
            if c > last and all(pc(c & b) <= 2 for b in F):
                rec(F + [c], c)
    rec([], -1)
    return classes

if __name__ == "__main__":
    if "--validate" in sys.argv:
        for n in (5, 6, 7):
            O = Orderly(n, 10 ** 6, "none"); O.Dmin = -10 ** 9; O.capd = 10 ** 9
            seen = []; O.check = lambda F, D: seen.append(tuple(sorted(F)))
            O.run()
            B = brute_classes(n, 100)
            assert len(set(seen)) == len(seen)
            print(f"validate n={n}: orderly classes={len(seen)} brute classes={len(B)} equal={set(seen) == B}")
        sys.exit()
    n, target, mode = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    O = Orderly(n, target, mode)
    print(f"n={n} target={target} mode={mode}: capd={O.capd} capl={O.capl} ell_max_global={O.ellmax_global} Dmin={O.Dmin}", flush=True)
    t0 = time.time(); res = O.run()
    print(f"nodes={O.nodes} families with D>=Dmin checked={O.families_checked} structures with count<={target}: {len(res)} "
          f"time={time.time()-t0:.1f}s")
    with open(f"orderly_n{n}_t{target}_{mode}.json", "w") as f:
        json.dump(res, f, indent=1)
