"""Independent exhaustive enumeration of 9-subsets of the 5x5 grid (audit of search-local).

Stage 1 (exact, integer/bitmask): for every 9-subset S of the 5x5 grid compute
    E(S)  = number of Euclidean circles through >= 3 points of S
    NB(S) = |B(S)| = number of Moebius blocks (lines or circles) through >= 3 points of S
    L(S)  = number of lines through >= 3 points of S
Stage 2: inversion-centre optimisation.  For a centre O not in S, the inverted set has
    circles = NB(S) - deg_S(O), deg_S(O) = number of blocks of S through O  (O = infinity: deg = L(S)).
Blocks of S through O pairwise share at most one point of S (they share O and at most one more),
so they form a partial linear space on 9 points with blocks of size >= 3: deg_S(O) <= 36/3 = 12.
Hence only subsets with NB(S) <= THRESH + 12 can reach a value <= THRESH; for those the exact
maximum of deg_S(O) is found by intersecting the blocks of S pairwise (every O of degree >= 2 is
such an intersection point), in floating point with residual-separation reporting, and the
near-record cases are re-verified exactly with sympy in grid5_exact_check.py.

Mode 'pure':      S = 9 grid points.
Mode 'with_inf':  S = 8 grid points + the point at infinity (then the centre O must be finite;
                  lines of the grid with >= 2 points of S8 are blocks of S).
"""
import sys, time, json
from itertools import combinations
from math import comb
import numpy as np
sys.path.insert(0, "/home/user/Solving-hard-math-problems-/experiments/506/audit-search-local")
from count_exact import block_key, on_block

K = 5
GRID = [(x, y) for x in range(K) for y in range(K)]
N = len(GRID)


def grid_blocks():
    blocks = {}
    for i, j, k in combinations(range(N), 3):
        key = block_key(GRID[i], GRID[j], GRID[k])
        if key not in blocks:
            blocks[key] = sum(1 << t for t, p in enumerate(GRID) if on_block(key, p))
    return blocks


def main(mode, thresh=31):
    t0 = time.time()
    blocks = grid_blocks()
    keys = list(blocks)
    masks = np.array([blocks[k] for k in keys], dtype=np.uint32)
    is_line = np.array([k[0] == 0 for k in keys])
    sizes = np.array([bin(int(m)).count("1") for m in masks])
    from collections import Counter
    print(f"[{mode}] grid blocks: {len(keys)} sizes={dict(sorted(Counter(sizes.tolist()).items()))} "
          f"lines={int(is_line.sum())} line sizes={dict(sorted(Counter(sizes[is_line].tolist()).items()))}")

    m = 9 if mode == "pure" else 8
    combs = np.array(list(combinations(range(N), m)), dtype=np.int8)
    S = np.zeros(len(combs), dtype=np.uint32)
    for c in range(m):
        S |= (np.uint32(1) << combs[:, c].astype(np.uint32))
    del combs
    print(f"[{mode}] subsets: {len(S)}  (expected C(25,{m}) = {comb(25, m)})  t={time.time()-t0:.1f}s")

    # correction tables: a block with k points of S accounts for C(k,3) triples but is one block
    corr3 = np.array([max(comb(k, 3) - 1, 0) for k in range(N + 1)], dtype=np.int32)   # blocks needing >=3
    corr2 = np.array([max(comb(k, 2) - 1, 0) for k in range(N + 1)], dtype=np.int32)   # lines needing >=2 (with_inf)
    hit3 = np.array([1 if k >= 3 else 0 for k in range(N + 1)], dtype=np.int32)
    hit2 = np.array([1 if k >= 2 else 0 for k in range(N + 1)], dtype=np.int32)

    E = np.zeros(len(S), dtype=np.int32)      # Euclidean circles
    NB = np.zeros(len(S), dtype=np.int32)     # |B(S)|
    L = np.zeros(len(S), dtype=np.int32)      # blocks through infinity
    for bi in range(len(keys)):
        k = np.bitwise_count(S & masks[bi]).astype(np.int32)
        if mode == "pure":
            if is_line[bi]:
                NB += hit3[k]; L += hit3[k]
            else:
                NB += hit3[k]; E += hit3[k]
        else:  # with_inf: infinity is in S, lines need >= 2 grid points, circles >= 3
            if is_line[bi]:
                NB += hit2[k]; L += hit2[k]
            else:
                NB += hit3[k]; E += hit3[k]
    if mode == "pure":
        # cross-check: E = C(9,3) - collinear triples - sum over circles (C(k,3)-1)
        E2 = np.full(len(S), comb(9, 3), dtype=np.int32)
        for bi in range(len(keys)):
            k = np.bitwise_count(S & masks[bi]).astype(np.int32)
            if is_line[bi]:
                E2 -= comb3(k)
            else:
                E2 -= corr3[k]
        assert np.array_equal(E, E2), "cross-check of Euclidean count failed"
    print(f"[{mode}] stage 1 done t={time.time()-t0:.1f}s")
    hist = lambda a: {int(v): int(c) for v, c in zip(*np.unique(a, return_counts=True))}
    out = {"mode": mode, "subsets": int(len(S)),
           "hist_euclid": hist(E), "hist_nb": hist(NB), "min_euclid": int(E.min()), "min_nb": int(NB.min())}
    if mode == "pure":
        print(f"[{mode}] Euclidean min = {E.min()}, count at min = {(E == E.min()).sum()}, "
              f"hist head = {sorted(out['hist_euclid'].items())[:10]}")
        for idx in np.nonzero(E == E.min())[0]:
            pts = [GRID[t] for t in range(N) if (int(S[idx]) >> t) & 1]
            print("   Euclidean minimiser:", pts, " NB=", NB[idx], " L=", L[idx])
            out.setdefault("euclid_minimisers", []).append(pts)
    print(f"[{mode}] min |B(S)| = {NB.min()}, hist_nb head = {sorted(out['hist_nb'].items())[:12]}")
    # lower bound for every centre of degree <= 12
    cand = np.nonzero(NB <= thresh + 12)[0]
    print(f"[{mode}] candidates with |B(S)| <= {thresh+12}: {len(cand)}")
    np.save(f"/home/user/Solving-hard-math-problems-/experiments/506/audit-search-local/cand_{mode}.npy",
            np.stack([S[cand].astype(np.int64), NB[cand], L[cand], E[cand]], axis=1))
    json.dump(out, open(f"/home/user/Solving-hard-math-problems-/experiments/506/audit-search-local/stage1_{mode}.json", "w"))
    json.dump({"keys": keys, "masks": [int(x) for x in masks]},
              open("/home/user/Solving-hard-math-problems-/experiments/506/audit-search-local/grid5_blocks.json", "w"))
    print(f"[{mode}] total t={time.time()-t0:.1f}s")


def comb3(k):
    return (k * (k - 1) * (k - 2)) // 6


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "pure")
