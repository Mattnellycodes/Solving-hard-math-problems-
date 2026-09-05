# audit-n10-enum / v2-cpsat — second independent audit of the n10-enum result (own code only)

Nothing here imports or reuses code from `n10-enum/`, from the earlier audit files in the parent directory, or from
`verify_independent/`; only the *data* files `n10-enum/orderly3/n10_sg_t32.json`, `n10_all4_sg_t32.json` are read for the
final isomorphism comparison.

## Method (deliberately different from the orderly generation used by n10-enum and by the first audit)
* `canon.py` — own canonical labelling of set systems (individualisation–refinement search tree, automorphism
  pruning, returns |Aut|).  Validated: the equivalence relation it induces coincides with brute-force S_n isomorphism
  on 600 random systems (n = 6, 7; 231 classes), |Aut| agrees, relabelling invariance on 200 random n = 10 systems,
  |Aut| of a 5-block / two disjoint 5-blocks / a 6-block = 14400 / 28800 / 17280 (`canon_validate.log`).
* `pipe.py` — the model.  Skeletons = families of blocks of size >= 5 pairwise sharing <= 2 points under the per-point
  SG cap, enumerated up to isomorphism by BFS + canon (**22 nonempty classes** for n = 10, same as the main agent's
  independent count).  `complete()` = OR-tools CP-SAT model with enumerate_all_solutions over the 4-blocks and the
  lines: every free triple in <= 1 rich block, 3-lines only on uncovered triples, lines pairwise <= 1 point,
  per-point derived cap C(n-1,2)-1, line cap C(n,2)-1, **plus the (valid) SG cap on the 10-point set Q_p**
  (images of P-p under inversion about p, plus p itself), D + ell >= C(n,3) - target.  `line_sets()` = own exact
  enumeration of all line sets of a rich-block structure.  `full_check()` = own hereditary checks on the 11 derived real
  point sets: A = hereditary Sylvester–Gallai (every non-collinear subset has an ordinary line), B = (8_3) on 8-subsets
  whose restricted lines are exactly 8 triples, K = Kelly–Moser (>= 3m/7 ordinary lines), C = cited o-table and
  orchard table t3.
* `planB.py n target i` — decomposition of skeleton i: pick a point q, O = its Aut(skeleton)-orbit, d = max 4-block
  degree on O; enumerate the derived structures at q (4-blocks through q) up to Aut(skeleton)_q-equivalence (canon with
  q individualised) and let CP-SAT complete each with no further 4-block through q and 4-degree <= d on O.  Needed
  because the direct CP-SAT enumeration of the single-5-block skeleton hit the 40-minute limit (UNKNOWN).
* `skel5.py` — special case of the same idea for the single-5-block skeleton (written first; agrees with planB).
* `all4.py` — the all-4-block case through a maximum-degree point: b4 >= 25, degrees <= 11, so the max degree d0 is
  10 or 11 and the derived structure at that point is a PSTS(9) with d0 triples; all PSTS(9) classes are enumerated
  (1,2,5,11,19,34,41,31,12,4,1 classes with 1..11 triples) and completed by CP-SAT.
* `theirs_check.py` — re-checks the 9 structures of their JSON with my line enumerator and my tiers.
* `mk_check.py` — (8_3): all 840 families of 8 pairwise <=1-sharing triples on 8 points form 1 class; with a frame on
  4 points containing no triple, the other points are forced and the closing condition is t^2 - t + 1 = 0: no real
  realisation with distinct points (independent of `verify_independent/mk_unrealisable2.py`).
* `compare.py` — isomorphism comparison of structure lists via my canon codes.

## Controls
* n = 8, target 18: big-block case empty; all-4-block case = 12-block class (count 17, 118 line sets) + 10-block class
  (count 18, 2 line sets, killed by (8_3)) — same as n10-enum.
* n = 9, target 25: big-block case = exactly their 2 classes ([5,5] sharing a point, count 24, 168 line sets, killed by
  hereditary SG; the 8-block, count 25, 105 line sets, passes everything); all-4-block case without the Q_p cap = exactly
  their 7 classes with line-set counts 3745/1/12/12/977/12/12 (with the cap: 5 of them); all die at (8_3).
* n = 9, target 27, single-5-block skeleton: `skel5.py` and `planB.py` give the same 9 classes (a direct CP-SAT
  enumeration of that skeleton was started as a third control and stopped after 55 min, unfinished).

## Result for n = 10, count <= 32 (D + ell >= 88), SG-only caps
* All 22 skeletons + the all-4-block case: every CP-SAT sub-problem finished with status OPTIMAL or INFEASIBLE.
* Big-block case: exactly **9 isomorphism classes**, pairwise isomorphic (my canon) to their 9 with identical counts and
  numbers of line sets: skeleton [6]: 2 classes (count 31, 28 line sets each); [5,5] sharing 2 points: 1 (count 32, 3);
  [5,5] disjoint: 3 (count 31; 30, 26, 31 line sets); [5,5,5] through a point: 1 (count 30, 342);
  [5,5,5,5]: 2 (count 32; 32 and 20). All-4-block case: exactly 1 class (27 four-blocks, count 32, 6 line sets).
* Own tier results on every line set: after hereditary SG + (8_3) the minimum count is 31; adding Kelly–Moser leaves only
  the two |Aut| = 40 classes of the disjoint-5-block skeleton with exactly 2 line sets each (count 32) and the two
  6-block classes with 2 line sets each (count 32); the orchard value t3(10) = 12 kills the latter.  The two surviving
  (F, L) pairs per class are exactly the line sets listed by n10-enum for its structures 5 and 6.

Logs: `run_n8_t18.log`, `run_n9_t25.log`, `all4_n8_t18.log`, `all4_n9_t25.log`, `all4_n9_t25_noQ.log`,
`skel5_n10_t32.log`, `planB_n10_t32_s*.log` (and `_noQ` variants), `all4_n10_t32.log`, `theirs_check.log`,
`mk_check.log`, `compare_*.log`, `canon_validate.log`; JSON outputs `out_*.json`.
