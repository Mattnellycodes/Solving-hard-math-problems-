# Erdős #506, n = 10: lower bound by CP-SAT (independent pipeline, `n10-cpsat/v2`)

Agent: n10-cpsat (v2), 2026-09-05.  Everything in this directory was written from scratch; no code from
`../` (the earlier n10-cpsat attempt), `../../n10-enum`, `../../theory` or `../../verify_independent` is imported.
The earlier attempt in `../` is kept untouched for reference; note that its CP-SAT enumerations stopped at their
time limits (status FEASIBLE), so they were not complete.

## 1. Statement and reduction

c(10) <= 32 would need an abstract Möbius structure (F, L) on P = {0..9} with D + l >= 88
(circles = C(10,3) - D - l, D = sum over rich blocks of C(|B|,3) - 1, l = number of lines).  By Lemmas A, B, C of
`../../theory/REPORT.md` the largest block has size 4, 5 or 6; the assigned task is the case "largest block 5 or 6"
(cases B and A below).  For completeness the pipeline also runs the case "largest block 4" (C) and the cases
"largest block 7, 8, 9" (D), so that the whole reduction of n = 10 is re-derived by the same model.

## 2. Facts used (all re-verified here unless marked CITED)

| tag | fact | where |
|---|---|---|
| T, L, PD | packing facts (triples in <= 1 rich block, lines pairwise <= 1 point, blocks through a pair partition the rest) | `model.py` |
| SG | derived Sylvester–Gallai: rich blocks through p cover <= C(9,2) - o(9) pairs; lines cover <= C(10,2) - o(10) | modes `sg` (o = 1, Sylvester–Gallai), `km` (o >= 3m/7, Kelly–Moser: o(9) >= 4, o(10) >= 5), `table` (o(9) = 6, o(10) = 5, Crowe–McKee, CITED) |
| HSG | hereditary SG on 7-/8-subsets of every derived point set and of P (excludes a derived Fano plane, F2) | `model.py`, `filters.py` |
| MK | 8 triples on 8 points pairwise sharing <= 1 point = Möbius–Kantor, unique (F3, 840 = 8!/48 labelled copies) and not realisable over R (residual equation t^2 - t + 1 = 0) | `facts.py` |
| O9 | 11 triples on 9 points pairwise sharing <= 1 point = AG(2,3) minus a line (all 10080 labelled systems checked) and contain an (8_3); hence t3(9) <= 10: at most 10 four-blocks through a point, at most 10 three-lines inside any 9 points | `facts.py` F4 |
| E11 | the same constraints for the derived structure at p of the 11-point Möbius set P ∪ {∞} (a rich line through p is ONE derived line, extended by ∞') | `model.py` |
| LPK | local packing on 9 points: with (d6, d5) derived 5-/4-lines at most 10 (MK), 8, 7, 4 / 6, 3 / 0 derived 3-lines; cuts 2 d4 + 3 d5 + 8 d6 <= 20, d5 + 2 d6 <= 3 | `facts.py` F1 + `model.py` |
| O10 | t3(10) <= 12: every 13-packing of triples on 10 points (leave graph = claw + 3 K_2, 144 labelled systems with the fixed leave, 2 isomorphism classes) is unrealisable over R (exact projective-frame construction; all solution components degenerate) | `facts_t3_10.py` (agrees with Burr–Grünbaum–Sloane 1974) |
| Miquel | for 8 distinct points of the Möbius plane labelled by a cube, the product of the six oriented edge cross-ratios is identically 1 (symbolic), so five concyclic faces force the sixth; used on P ∪ {∞} | `facts.py` F5, `filters.py` |
| (not used) | the naive bundle theorem is FALSE in degenerate position (the antipodal 10-set violates it) | `filters.py` docstring |

Validation of the model: n = 8 gives optimum D + l = 39, i.e. c(8) >= 17 with the cube structure (`validate.log`),
matching the established value.

## 3. Pipeline

1. `skeletons.py`: isomorphism classes of families of big blocks (sizes 5, 6; <= 2 six-blocks and <= 3 five-blocks
   through a point by F1): 6 skeletons containing a 6-block (case A), 12 containing a 5-block and no 6-block (case B);
   plus the single skeletons for cases C (4-block) and D (7-, 8-, 9-block).
2. `run_skeletons.py`: CP-SAT feasibility of D + l >= 88 with the skeleton's big blocks fixed exactly (other big
   blocks forbidden), symmetry broken by non-increasing degree keys inside the cells of the skeleton.
3. `enum_skeleton.py`: for feasible skeletons, complete enumeration of the 4-block families by CP-SAT with
   *orbit nogoods* (after each solution its whole orbit under Aut(skeleton) is excluded by exact-copy clauses;
   termination = INFEASIBLE = completeness), then exhaustive enumeration of all line sets per 4-block family,
   reduction of (F, L) to isomorphism classes.
4. `postprocess.py`: hereditary SG / hereditary orchard / Miquel filters (`filters.py`), numerical realisation
   (`realise.py`, random-start least squares with structure verification).
5. Exact analyses of the survivors: `analyse_six.py` (6-block: pencil/Frégier involution conditions, Groebner),
   `analyse_two5.py` (two disjoint 5-blocks: Möbius normal forms, torsion lattice of the angle system, Groebner
   for the point at infinity).

## 4. Results

(filled in below as the runs complete)

### 4.0 A hand proof for the single-5-block skeleton (case B, k = 1), using O10

Let B5 = {0..4} be the only block of size >= 5.  Then D = 9 + 3 b4 and, with l <= 14 (SG, o(10) >= 1),
D + l >= 88 forces b4 >= 22.  At p ∈ B5 the derived structure contains the 4-line B5 - p, so d4(p) <= 8 (F1); at
p ∉ B5, d4(p) <= 10 (O9).  Hence 4 b4 = sum_p d4(p) <= 5·8 + 5·10 = 90, so b4 = 22, every cap is attained
(d4 = 8 on B5, 10 off B5), D = 75 and l >= 13.
Let a_p, b_p be the numbers of 3-lines and 4-lines through p (a 4-line is a 4-block).  In the derived 10-point real
set Q_p at p (E11) the 3-point lines are the d4(p) - b_p derived 3-lines not extended by ∞' plus the a_p lines
(S - p) ∪ {∞'}; by O10 (t3(10) <= 12): a_p - b_p <= 2 (p ∉ B5), a_p - b_p <= 4 (p ∈ B5).  Pair counting in Q_p
(SG with o = 1: <= 44 covered pairs; a rich 4-line through p covers 6 pairs, a derived 3-line 3, the derived 4-line
B5 - p 6 or, if B5 is a line, 10): a_p + b_p <= 4 at every point.  For p ∉ B5 this gives a_p + b_p <= 2 + 2 b_p,
i.e. a point off B5 lies on at most 2 lines unless it is on a 4-line.  Now
sum_p (a_p + b_p + [p ∈ B5 and B5 a line]) = sum_S |S| = 3 l + #4-lines + 2·#5-lines >= 39 + #4-lines + 2·#5-lines,
while the left side is <= 5·4 + 5·4 = 40.  Hence #4-lines + 2·#5-lines <= 1.  If there is no 4-line, the five points
off B5 contribute <= 2 each and the total is <= 30 < 39, contradiction.  If there is exactly one 4-line (and no
5-line) then the total must be exactly 40 = 39 + 1, so every point has exactly 4 lines, so every point off B5 lies
on a 4-line — but the unique 4-line has only 4 points.  Contradiction.  So no structure with largest block 5 and a
single 5-block reaches D + l >= 88.  (The CP-SAT run proves the same without O10 if it terminates; see 4.1.)

### Implementation notes (bugs found and fixed during development)
* First version of the E11 constraints double-counted a rich line through p (as B - p and as (B - p) ∪ {∞'}),
  which made the constraints too strict (invalid).  Fixed before the runs reported below (the first, invalid
  run was killed; its log is `skel_sg_buggyE11_killed.log`).  The same double counting was then found and fixed
  in `filters.derived11`; a consistency test (the CP-SAT solution for skeleton 0 must pass every filter that the
  model encodes) is what caught it.
* The naive bundle filter was removed after it flagged the (realisable) antipodal configuration.
* Tests: `analyse_two5.py` on the concentric-pentagon block family with a realisable 2-line set returns a
  non-trivial ideal (finite O) and (1) for O = ∞, as it must; with a full 10-line set it returns (1) everywhere.
  `analyse_six.py` on the skeleton-0 example structure returns (1) (6 involution conditions).
* Robustness: the case-A skeletons 1–5 are infeasible also WITHOUT the E11 constraints (`skel_sg_noe11_A.log`),
  and z3 (independent pseudo-Boolean encoding, `z3_check.py`) confirms infeasibility of skeleton 17 in 4 s
  (`z3_easy.log` for the others).
