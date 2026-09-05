# search-local — file index (Erdős #506, local-search construction agent, 2026-09-05)

Code (all independent of the other agents' code except the exact certificate cross-check):
- `engine.py`      — universe/block structure (exact integer block finder for rational universes,
                     float finder with neighbour-merged quantisation otherwise), numba incremental
                     counter (add/remove a point, |B(S)| and per-point degrees), simulated annealing
                     with tabu list and degree-biased candidate choice (`k_sa`), lexicographic
                     exhaustive DFS (`k_exhaustive`), `true_count` (best inversion centre over all
                     pairwise block intersections).
- `universes.py`   — builders: grids, grid + inversions in lattice points, triangular lattice,
                     regular polygons at several radii, orthocentric closures (rational), icosahedral /
                     octahedral sphere point sets (stereographic projection), ellipse (eccentric-angle
                     orchard), rectangular hyperbola, regular-polygon diagonal closures,
                     arrangement points of multiplicity >= 3 (`multiplicity_points`).
- `run_sa.py`      — SA driver over the universe catalogue (`--closure` adds the highest-multiplicity
                     arrangement points as selectable points); appends to `runs/sa_results.jsonl`,
                     records with count <= f(n) to `runs/sa_records.jsonl`.
- `exhaustive_grid.py` — exhaustive n-subsets of the k x k grid with the inversion centre ranging over
                     grid points, infinity and all multiplicity>=3 arrangement points (exact treatment
                     of every centre of degree >= 3; |B(S)| histogram covers degree <= 2 centres).
- `verify_exact.py` — exact certification (Fractions + cross-check with ../circles_exact.py for
                     rational sets; sympy-certified coincidences + 60-digit separation otherwise);
                     identifies the inversion centre exactly as an intersection of two blocks.
- `make_records.py`, `antipodal_rational.py`, `summarise.py`, `test_engine.py`.

Runs: `runs/batch*.sh` + logs, `runs/exhaustive_grid5_n9.{log,json}`, `runs/sa_results.jsonl`,
`runs/sa_records.jsonl`, `runs/records.json`, `runs/antipodal_rational.json`.
