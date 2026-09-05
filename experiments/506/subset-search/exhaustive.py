"""Exhaustive enumeration of all n-subsets of a k x k integer grid (Erdos #506).
For every subset we record (a) the Euclidean circle count of the lattice configuration itself
(O = infinity) and (b) the Moebius-optimised count with the inversion centre restricted to grid
points (plus infinity); infinity is NOT a member of the subsets here (pure lattice configurations),
except in the '--with-inf' variant.  Degenerate subsets (all on one block) are skipped."""
import numpy as np, sys, time, json, itertools
from numba import njit
import universe
from search import _add, _remove, BIG

@njit(cache=True)
def enum(n, N1, allow_inf, bptr, bmem, pptr, pblk):
    M = bptr.shape[0] - 1
    inS = np.zeros(N1, np.bool_)
    cnt = np.zeros(M, np.int32)
    deg = np.zeros(N1, np.int32)
    S = np.empty(n, np.int32)
    top = N1 if allow_inf else N1 - 1
    best_e = 10**9; best_m = 10**9
    bestS_e = np.zeros(n, np.int32); bestS_m = np.zeros(n, np.int32); bestO_m = -1
    hist_e = np.zeros(400, np.int64); hist_m = np.zeros(400, np.int64)
    nb = 0; full = 0
    depth = 0
    S[0] = -1
    total = 0
    # iterative DFS over combinations in lexicographic order
    while depth >= 0:
        S[depth] += 1
        if S[depth] > top - (n - depth):
            depth -= 1
            if depth >= 0:
                a, b = _remove(S[depth], n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
            continue
        a, b = _add(S[depth], n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
        if depth == n - 1:
            total += 1
            if full == 0:
                # Euclidean count: nblocks - deg[inf]
                e = nb - deg[N1 - 1]
                md = 0; mo = -1
                for p in range(N1):
                    if not inS[p] and deg[p] > md:
                        md = deg[p]; mo = p
                m = nb - md
                if e < 400: hist_e[e] += 1
                if m < 400: hist_m[m] += 1
                if e < best_e:
                    best_e = e; bestS_e[:] = S
                if m < best_m:
                    best_m = m; bestS_m[:] = S; bestO_m = mo
            a, b = _remove(S[depth], n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
        else:
            depth += 1
            S[depth] = S[depth - 1]
    return total, best_e, bestS_e, best_m, bestS_m, bestO_m, hist_e, hist_m

def run(k, n, allow_inf=False):
    U = universe.grid(k)
    t = time.time()
    total, be, Se, bm, Sm, Om, he, hm = enum(n, U.N + 1, allow_inf, U.bptr, U.bmem, U.pptr, U.pblk)
    dt = time.time() - t
    from search import formula
    r = dict(grid=f"{k}x{k}", n=n, allow_inf=allow_inf, subsets=int(total), seconds=round(dt, 1), formula=formula(n),
             min_euclid=int(be), argmin_euclid=[U.exact_str(i) for i in Se],
             min_moebius=int(bm), argmin_moebius=[U.exact_str(i) for i in Sm], O=U.exact_str(Om) if Om >= 0 else None,
             hist_euclid={i: int(he[i]) for i in range(400) if he[i]}, hist_moebius={i: int(hm[i]) for i in range(400) if hm[i]})
    print(json.dumps(r))
    sys.stdout.flush()
    return r

if __name__ == "__main__":
    jobs = [(4, 6), (4, 7), (4, 8), (4, 9), (4, 10), (5, 6), (5, 7), (5, 8), (5, 9), (5, 10)]
    if len(sys.argv) > 1:
        jobs = [tuple(map(int, a.split(','))) for a in sys.argv[1:]]
    res = []
    for k, n in jobs:
        res.append(run(k, n, False))
        if k == 4 or n <= 8:
            res.append(run(k, n, True))
    json.dump(res, open(f"exhaustive_{'_'.join(f'{k}x{k}n{n}' for k,n in jobs)}.json", "w"), indent=1)
