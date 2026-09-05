"""Exhaustive search over all n-subsets of the k x k grid (Erdős #506, search-local agent).

Two problems:
 (a) pure lattice subsets S ⊆ grid (|S| = n): Euclidean circle count (O = infinity) and the
     Möbius-optimised count min over ALL inversion centres O ∉ S.  O ranges over the grid points,
     infinity and every point where >= 3 blocks of the grid universe meet (these are all the points
     that can have deg_S(O) >= 3); a centre of degree <= 2 can only give |B(S)| - 2, which is checked
     against the histogram of |B(S)|.
 (b) S ⊆ grid ∪ {infinity} (n points) — i.e. n-1 grid points plus infinity, inverted about a grid
     point / arrangement point.
"""
import sys, time, json
import numpy as np
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/search-local')
from engine import *
import universes as UV

k = int(sys.argv[1]); n = int(sys.argv[2])
t0 = time.time()
G = UV.grid(k)
print(G.summary(), flush=True)
ext = UV.multiplicity_points(G, min_mult=3, max_pairs=10**9)
print(f"arrangement points of multiplicity>=3 (not grid points): {len(ext)}  "
      f"(max mult {max(e[2] for e in ext) if ext else 0})  t={time.time()-t0:.1f}s", flush=True)
E = np.array([[e[0], e[1]] for e in ext], dtype=np.float64) if ext else None
U = UV.grid(k, name=f"grid{k}x{k}+arr")
if E is not None:
    U.add_extra_O(E)
    U._csr()
print(U.summary(), "T =", U.T, flush=True)
out = {}
for mode, M, einf in (("pure", U.N, U.N), ("with_inf", U.M, U.N)):
    t = time.time()
    total, best, hist, hist_e, hist_nb, good = k_exhaustive(n, M, U.T, U.bptr, U.bmem, U.pptr, U.pblk, best_thresh := formula(n), einf)
    dt = time.time() - t
    h = {int(i): int(hist[i]) for i in range(400) if hist[i]}
    he = {int(i): int(hist_e[i]) for i in range(400) if hist_e[i]}
    hnb = {int(i): int(hist_nb[i]) for i in range(400) if hist_nb[i]}
    min_nb = min(hnb) if hnb else None
    print(f"[{mode}] n={n} subsets={total} best_moebius={best} best_euclid={min(he) if he else None} "
          f"min|B(S)|={min_nb} (so deg<=2 centres give >= {min_nb-2}) time={dt:.1f}s", flush=True)
    print("   hist moebius (first 8):", sorted(h.items())[:8])
    print("   hist euclid  (first 8):", sorted(he.items())[:8])
    # verify the good subsets with the true count (all O) and describe
    recs = []
    for S in good[:200]:
        tc, O, d, nb, lines = true_count(U, S)
        recs.append(dict(S=[U.label(int(i)) for i in S], proxy=int(k_eval(np.array(S), n, U.T, U.bptr, U.bmem, U.pptr, U.pblk)[0]),
                         true=int(tc), centre=(None if O is None else (O if isinstance(O, str) else [float(O[0]), float(O[1])])),
                         structure=describe(U, S)))
    if recs:
        print("   example:", recs[0])
    out[mode] = dict(n=n, k=k, subsets=int(total), best_moebius=int(best), best_euclid=(min(he) if he else None),
                     min_nb=min_nb, hist_moebius=h, hist_euclid=he, hist_nb=hnb, seconds=dt,
                     n_good=int(len(good)), good_examples=recs[:50])
json.dump(out, open(f"runs/exhaustive_grid{k}_n{n}.json", "w"), indent=1)
print("total time", time.time() - t0)
