"""Simulated annealing over n-subsets of a universe (see universe.py), minimising
   value(S) = |B(S)| - max_{O in U' \ S} deg_S(O)   (+ penalty if S lies on one block),
i.e. the number of circles of the best inversion of S with the inversion centre restricted to
universe points (or infinity).  Incremental counter in numba."""
import numpy as np, math, time, json, sys, os
from numba import njit

BIG = 10_000

@njit(cache=True)
def _add(p, n, inS, cnt, deg, bptr, bmem, pptr, pblk):
    """add point p to S; returns (delta_nblocks, delta_full) where full = #blocks with cnt == n"""
    dnb = 0; dfull = 0
    inS[p] = True
    for t in range(pptr[p], pptr[p + 1]):
        b = pblk[t]
        c = cnt[b] + 1
        cnt[b] = c
        if c == 3:
            dnb += 1
            for u in range(bptr[b], bptr[b + 1]):
                deg[bmem[u]] += 1
        if c == n:
            dfull += 1
    return dnb, dfull

@njit(cache=True)
def _remove(p, n, inS, cnt, deg, bptr, bmem, pptr, pblk):
    dnb = 0; dfull = 0
    inS[p] = False
    for t in range(pptr[p], pptr[p + 1]):
        b = pblk[t]
        c = cnt[b] - 1
        if cnt[b] == n:
            dfull -= 1
        cnt[b] = c
        if c == 2:
            dnb -= 1
            for u in range(bptr[b], bptr[b + 1]):
                deg[bmem[u]] -= 1
    return dnb, dfull

@njit(cache=True)
def _value(nb, full, inS, deg):
    md = 0; mo = -1
    for p in range(inS.shape[0]):
        if not inS[p] and deg[p] > md:
            md = deg[p]; mo = p
    v = nb - md
    if full > 0:
        v += BIG
    return v, mo

@njit(cache=True)
def sa_run(n, N1, bptr, bmem, pptr, pblk, iters, T0, T1, seed, init, tabu_len):
    np.random.seed(seed)
    M = bptr.shape[0] - 1
    inS = np.zeros(N1, np.bool_)
    cnt = np.zeros(M, np.int32)
    deg = np.zeros(N1, np.int32)
    S = np.empty(n, np.int32)
    nb = 0; full = 0
    if init[0] >= 0:
        for i in range(n):
            S[i] = init[i]
    else:
        perm = np.random.permutation(N1)
        for i in range(n):
            S[i] = perm[i]
    for i in range(n):
        a, b = _add(S[i], n, inS, cnt, deg, bptr, bmem, pptr, pblk)
        nb += a; full += b
    cur, curO = _value(nb, full, inS, deg)
    best = cur; bestS = S.copy(); bestO = curO
    recent = np.full(tabu_len, -1, np.int32); rp = 0
    lT0 = math.log(T0); lT1 = math.log(T1)
    for it in range(iters):
        T = math.exp(lT0 + (lT1 - lT0) * it / iters)
        i = np.random.randint(n)
        p = S[i]
        q = np.random.randint(N1)
        if inS[q]:
            continue
        # tabu: do not re-add a recently removed point
        tab = False
        for t in range(tabu_len):
            if recent[t] == q:
                tab = True; break
        if tab:
            continue
        a, b = _remove(p, n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
        a, b = _add(q, n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
        new, newO = _value(nb, full, inS, deg)
        d = new - cur
        if d <= 0 or np.random.random() < math.exp(-d / T):
            S[i] = q; cur = new; curO = newO
            if tabu_len > 0:
                recent[rp] = p; rp = (rp + 1) % tabu_len
            if cur < best:
                best = cur; bestS[:] = S; bestO = curO
        else:
            a, b = _remove(q, n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
            a, b = _add(p, n, inS, cnt, deg, bptr, bmem, pptr, pblk); nb += a; full += b
    return best, bestS, bestO

@njit(cache=True)
def eval_subset(S, N1, bptr, bmem, pptr, pblk):
    n = S.shape[0]
    M = bptr.shape[0] - 1
    inS = np.zeros(N1, np.bool_)
    cnt = np.zeros(M, np.int32)
    deg = np.zeros(N1, np.int32)
    nb = 0; full = 0
    for i in range(n):
        a, b = _add(S[i], n, inS, cnt, deg, bptr, bmem, pptr, pblk)
        nb += a; full += b
    v, o = _value(nb, full, inS, deg)
    return v, o, nb

def formula(n):
    return (n - 1) * (n - 2) // 2 + 1 - (n - 1) // 2

def run_universe(U, ns, runs=4, iters=400_000, T0=2.0, T1=0.05, seed=0, tabu_len=3, log=print):
    """returns dict n -> list of (value, sorted subset, O) best results (distinct subsets)."""
    N1 = U.N + 1
    out = {}
    for n in ns:
        if n >= N1: continue
        res = {}
        for r in range(runs):
            init = np.full(n, -1, np.int32)
            best, bestS, bestO = sa_run(n, N1, U.bptr, U.bmem, U.pptr, U.pblk, iters, T0, T1, seed * 1000 + r * 17 + n, init, tabu_len)
            key = tuple(sorted(int(x) for x in bestS))
            if best < BIG:
                res[key] = (int(best), int(bestO))
        if not res:
            continue
        items = sorted(((v, k, o) for k, (v, o) in res.items()))
        out[n] = items
        bv = items[0][0]
        flag = " <= formula" if bv <= formula(n) else ""
        log(f"  {U.name} n={n}: best {bv} (formula {formula(n)}) O={items[0][2]}{flag}  distinct={len(items)}")
    return out

if __name__ == "__main__":
    import universe
    U = universe.grid(5)
    t = time.time()
    r = run_universe(U, [8], runs=2, iters=50_000)
    print(time.time() - t)
    t = time.time()
    r = run_universe(U, range(6, 16), runs=3, iters=300_000)
    print(time.time() - t)
