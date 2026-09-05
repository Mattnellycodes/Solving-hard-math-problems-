# Erdős #506, n = 10: c(10) >= 33 (hence c(10) = 33 = f(10)) — n10-continue agent, 2026-09-05

Directory `experiments/506/n10-continue/`.  All code here was written from scratch (nothing imported from
`n10-cpsat/v2`, `n10-enum`, `theory`, `verify_independent`); their results are used only as cross-checks.

## 0. Result

**Theorem.** Every set of 10 points in the plane, not all on one line or one circle, determines at least 33
circles through >= 3 of its points.  With the antipodal configuration (8 concyclic points in 4 antipodal pairs +
the centre + one more point on a diameter, 33 circles, `../theory/certify.py`) this gives **c(10) = 33 = f(10)**.

**Assumptions (all theorems; no cited tables are used in the main line):**
* (SG) Sylvester–Gallai: a finite non-collinear real point set has an ordinary line (o(m) >= 1);
* (8_3) the Möbius–Kantor configuration (8 points, 8 collinear triples pairwise sharing <= 1 point) has no
  realisation by 8 distinct real points — proved from scratch in this session (`../n10-cpsat/v2/facts.py` F3,
  `../verify_independent/mk_unrealisable2.py`, adjudicated in `../skeptics/ADJUDICATION.md`);
* (STS9) every 11 triples on 9 points pairwise sharing <= 1 point contain an (8_3) — re-verified here
  (`t3_9_check.py`: 10080 labelled 11-packings, all contain an (8_3)); with (8_3) this gives t3(9) <= 10;
* (MIQ) Miquel's theorem in the form "8 distinct points of the Möbius plane labelled by a cube never have exactly
  five concyclic faces" — the underlying cross-ratio identity (product of the six face cross-ratios = 1) is
  re-verified symbolically here (`miquel_identity.py`);
* elementary packing/parity facts proved in the text (t3(10) <= 13 by parity; per-point caps);
* exact polynomial algebra (Gröbner bases, sympy 1.14) with explicit certificates printed in the logs;
* the correctness of the enumerations (CP-SAT 9.15 with all-solution enumeration; own DFS), cross-checked as
  described in §6.
The cited values o(9) = 6, o(10) = 5 (Crowe–McKee) and t3(10) = 12 (Burr–Grünbaum–Sloane) are NOT needed;
they are only reported as an additional, independent kill of the same structures (§5).

## 1. Setting and reduction (Möbius reformulation, as in `../theory/REPORT.md` §1)

P = {0..9}.  Blocks = circles-or-lines through >= 3 points; every triple lies in exactly one block; two blocks
share <= 2 points; lines (blocks through ∞) pairwise share <= 1 point.  circles(P) = 120 − D − l with
D = Σ_{|B|>=4} (C(|B|,3) − 1), l = #lines (size >= 3).  circles(P) <= 32  ⟺  D + l >= 88.
Inverting P ∪ {∞} about p ∈ P gives a real 10-point set Q_p = (P − p) ∪ {∞'} (not collinear, since P is not
collinear) whose >= 3-point lines are B − p for the rich blocks B ∋ p (extended by ∞' if B is a line) and
(S − p) ∪ {∞'} for the 3-lines S ∋ p; dropping ∞' gives the 9-point derived set.  By SG the derived 9-set's
lines cover <= 35 pairs, Q_p's lines <= 44 pairs, the lines of P <= 44 pairs (so l <= 14).

### 1.1 Largest block >= 7 is excluded (`lemmasABC.py`, re-derived)
* m = 9 (Lemma A): the 36 triples {x, a, b} (a, b ∈ B) lie in distinct blocks (a common block would share 3 points
  with B); if B is a line no other line exists (it would share 2 points with B), else lines through x cut B in
  disjoint pairs, l <= 4: circles >= 37 − 4 = 33.
* m = 8 (Lemma B): with x, y ∉ B and t <= 4 four-blocks {x, y, a, b}: |blocks| = 65 − 3t >= 53, l <= 14:
  circles >= 39.
* m = 7 (Lemma C, chunk lemma): N2 >= r·C(m,2) − m·r(r−1)/4 = 63 − 10.5 → >= 53 blocks meeting B in two points:
  circles >= 54 − 14 = 40.  (Proof of the chunk lemma checked line by line: for fixed a ∈ B the chunks of
  different b share <= 1 point, so Σ_b Σ_chunks C(|S|,2) <= C(r,2), and r − g_ab <= Σ_chunks C(|S|,2).)
So D + l >= 88 forces a largest block of size 4, 5 or 6 (no rich block at all gives circles >= 120 − 14).

### 1.2 Largest block 4 (all blocks of size <= 4) — hand proof, theorems only
b4 = #4-blocks.  Each point lies in <= t3(9) = 10 four-blocks (their derived 3-point lines form a packing on 9
points; (STS9)+(8_3)), so b4 <= 25; D = 3 b4 and l <= 14 give b4 >= 25.  Hence b4 = 25, every point in exactly
10 four-blocks, D = 75, l >= 13.  Let a_p, b_p be the numbers of 3-lines and 4-lines through p.  In Q_p the
lines are 10 − b_p three-point lines (non-line 4-blocks), b_p four-point lines, a_p three-point lines through ∞';
SG on Q_p: 3(10 − b_p) + 6 b_p + 3 a_p <= 44, i.e. a_p + b_p <= 4.  Summing, 3 l3 + 4 l4 = Σ_p (a_p + b_p) <= 40
with l3 + l4 >= 13, so l4 <= 1.  Parity: 14 lines with exactly 3 points on 10 real points are impossible (they would
cover 42 pairs, the 3 uncovered pairs must give every point an odd leave-degree, 10 odd degrees need >= 10
endpoint incidences > 6).  If l4 = 0: Σ a_p = 3 l >= 39 forces a point with a_p = 4, whose Q_p has 10 + 4 = 14
three-point lines — contradiction.  If l4 = 1: Σ = 40 forces a_p + b_p = 4 at every point, and a point off the
4-line has b_p = 0, a_p = 4, again 14 three-point lines in Q_p — contradiction.  (Note: the earlier reports
argued "l <= 12 by t3(10) = 12", which overlooks l = 12 three-lines + one 4-line; the argument above closes
that gap and needs no cited value.)

### 1.3 Skeletons (largest block 5 or 6) — own enumeration (`skeletons_own.py`)
All families of blocks of size 5..9 on 10 points, pairwise sharing <= 2 points, under the per-point SG cap
Σ_{B∋p} C(|B|−1, 2) <= 35: **22 isomorphism classes** (BFS with invariant + VF2), identical to
`../verify_independent/n10/skeletons_check.py`; they consist of the 18 skeletons of `../n10-cpsat/v2/skeletons.json`
(largest block 5 or 6; the v2 list is complete) plus [7], [8], [9], [5,7], which contain a block of size >= 7
and are excluded by §1.1.  (Note: the v2 list did not contain [5,7]; harmless, as [5,7] is excluded by Lemma C.)

## 2. Method for the 18 skeletons (`enum4.py`)

For each skeleton (the exact family of blocks of size >= 5):
* **Stage 1** (CP-SAT, complete enumeration of all solutions): 4-block families F4 with (T) every triple in <= 1
  block; (SG) per-point cap 35; (MK) for every p and every 8-subset E of the other points, at most 7 blocks
  B ∋ p with |(B − p) ∩ E| = 3 (the (8_3) lemma; the restriction of a derived 4- or 5-line counts); (t3(9))
  <= 10 four-blocks through p; (LB) D_skel + 3 b4 >= 88 − 14.  Symmetry breaking by lex-leader constraints
  x <=_lex g(x) for a generating set of Aut(skeleton) plus all its transpositions (valid: the lex-smallest member of
  each orbit satisfies them; unit-tested).  Solutions reduced to isomorphism classes (invariant + VF2).
* **Stage 2** (own DFS, exhaustive): all line sets L with |L| >= 88 − D(F): rich lines are blocks of F, 3-lines are
  uncovered triples, lines pairwise share <= 1 point, Σ C(|S|,2) <= 44.  (F, L) reduced to isomorphism classes.
* **Filters** (`postfilter.py`, `common.py`): for the 21 real linear spaces forced by (F, L) — the derived 9-sets,
  the 10-sets Q_p, and (P, L) — hereditary SG with the correct hypothesis (a subset contained in one line is
  skipped), orchard caps t3(7..9) = 6, 7, 10 and t3(10) <= 13, and Miquel closure on P ∪ {∞}.
* **Exact analysis** of what survives: §4.

## 3. Results per skeleton (`runs/summary_table.md`, `runs/enum_*.log`, `runs/post_all.log`)

| own | v2 | big blocks | D_skel | b4_min | |Aut| | stage 1 | labelled | F-classes | (F,L) classes | theorem-only kills |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 6 | [5] | 9 | 22 | 14400 | INFEASIBLE | 0 | 0 | 0 | — |
| 1 | 0 | [6] | 19 | 19 | 17280 | complete | 30 | 8 | 8 | Miquel x8 (+ exact: §4.1) |
| 5 | 7 | [5,5] (2 common) | 18 | 19 | 288 | complete | 30 | 6 | 0 | no line set reaches l >= 88 − D |
| 6 | 8 | [5,5] (1 common) | 18 | 19 | 1152 | INFEASIBLE | 0 | 0 | 0 | — |
| 7 | 9 | [5,5] disjoint | 18 | 19 | 28800 | complete | 3054 | 59 | 16 | Miquel x14, orchard x1, 1 survivor -> exact (§4.2); exact kills all 16 |
| 8 | 1 | [5,6] | 28 | 16 | 288 | INFEASIBLE | 0 | 0 | 0 | — |
| 9 | 2 | [5,6] | 28 | 16 | 2880 | INFEASIBLE | 0 | 0 | 0 | — |
| 11 | 3 | [6,6] | 38 | 12 | 2304 | complete | 1 | 1 | 0 | no line set |
| 12 | 10 | [5,5,5] | 27 | 16 | 48 | complete | 8 | 2 | 0 | no line set |
| 13 | 11 | [5,5,5] | 27 | 16 | 48 | complete | 34 | 10 | 0 | no line set |
| 14 | 12 | [5,5,5] | 27 | 16 | 32 | INFEASIBLE | 0 | 0 | 0 | — |
| 15 | 4 | [5,5,6] | 37 | 13 | 32 | INFEASIBLE | 0 | 0 | 0 | — |
| 16 | 13 | [5,5,5,5] | 36 | 13 | 8 | complete | 1 | 1 | 0 | no line set |
| 17 | 14 | [5,5,5,5] | 36 | 13 | 8 | INFEASIBLE | 0 | 0 | 0 | — |
| 18 | 5 | [5,5,5,6] | 46 | 10 | 48 | INFEASIBLE | 0 | 0 | 0 | — |
| 19 | 15 | [5,5,5,5] | 36 | 13 | 128 | complete | 18 | 9 | 13 | hereditary SG (derived Fano) x13 (also orchard, Miquel) |
| 20 | 16 | [5^5] | 45 | 10 | 10 | INFEASIBLE | 0 | 0 | 0 | — |
| 21 | 17 | [5^6] | 54 | 7 | 60 | INFEASIBLE | 0 | 0 | 0 | — |

Total: 37 (F, L) classes with D + l >= 88 over all skeletons (`runs/post_all.json`); 36 are killed by the
theorem-only filters; the single survivor (structure 12 of `post_all`, = structure 4 of `post_7`) is killed by
the exact analysis of §4.2.  All 18 skeletons are therefore excluded, with §1.1–1.3 this proves the theorem.

### 3.1 Skeleton [5] (single 5-block; v2 skeleton 6, "UNKNOWN after 2400 s" there)
Stage 1 is INFEASIBLE in 4 s (41 s without symmetry breaking, `runs/enum_0_nolex.log`; also without the explicit
t3(9) cap, `runs/enum_0_nocap.log`).  What kills it: b4 >= 22 (from l <= 14), and no family of 22 four-blocks
exists once the (8_3) constraint is imposed at the points of the 5-block *including the restriction of the derived
4-line B5 − p to 8-subsets*.  **The hand proof in `../n10-cpsat/v2/REPORT.md` §4.0 is wrong**: it claims that
4 b4 = Σ_p d4(p) <= 5·8 + 5·10 = 90 forces "b4 = 22 and every cap attained", but 4·22 = 88 < 90 leaves slack 2
(a CP-SAT solution with the 5-block, 22 four-blocks pairwise sharing <= 2 points, d4 = [8,7,8,8,8,10,9,10,10,10],
SG cap and d4 <= 10 off the block exists: `skel6_inspect.py`).  Its conclusion is nevertheless true, by the
computation above.  A weaker but sufficient hand argument is not attempted here.

## 4. Exact analyses of the survivors

### 4.1 Single 6-block (skeleton [6]; the 8 F-classes have b4 = 19, degrees 7 on the 6-block and 10 off it)
Only F-classes 0 and 6 admit line sets (28 labelled each; 8 (F, L) classes, l = 12 or 13, counts 32/31).
All 8 (F, L) classes violate Miquel closure (`runs/post_all.log`).  Independently (`analyse_six_own.py`,
`involution_certificates.py`): for r, r' ∉ A the circles through r, r' cut the circle A in pairs of one projective
involution of A (pencil/Frégier argument re-derived in the docstring; the centre is not on A because three blocks
through r, r' cut A in disjoint pairs), so each pair of R gives a 3x3 determinant condition on the parameters
(0, 1, ∞, t3, t4, t5) of the six points of A.  Saturated by distinctness:
F-class 0: Gröbner basis (1) (unrealisable over C);  F-class 6: lex basis [t3 − t5 − 1, 2 t4 − t5 − 1, t5² + 1]
(no real point).  The other six F-classes (no line sets anyway): classes 4 → (1); 1, 2, 3, 7 → only complex
points (t5² + 1 or 2 t5² − 2 t5 + 1); class 5 → a real 1-parameter family (the anharmonic 6-tuple) but no line set.

### 4.2 Two disjoint 5-blocks (skeleton [5,5] disjoint; `analyse_two5_own.py`, `two5_families.py`)
Every 4-block is {s, s', r, r'}.  With S on circle A and R on circle B (A ≠ B, no point of P on both) a Möbius
normal form gives (i) concentric circles |z| = 1, |z| = ρ: concyclic ⟺ α_s + α_s' ≡ β_r + β_r' (mod 2π)
(common perpendicular bisector of the two chords); (ii) parallel lines: x_s + x_s' = t_r + t_r'; (iii) two lines
through 0: x_s x_s' = t_r t_r'.  So the 4-blocks give an integer linear system M θ = 0.  For all 59 F-classes the
real kernel of M is spanned by (1,…,1), which excludes (ii), (iii) (points of S would coincide / have equal
modulus), and in (i) every solution is a rotation of a torsion solution (Smith normal form, own implementation,
tested against brute force).  56 F-classes have torsion exponent 2 and no labelling with 5 distinct points on a
circle: their 4-block families are unrealisable whatever the lines.  The 3 remaining F-classes (19, 20, 54) have
e = 10 and 8 labellings each: they are the concentric regular pentagons (angles 2πk/10; F-class 19 = the
well-known 20-circle family, the only one with line sets: 31 labelled, 5 (F, L) classes).  For every (F, L) class
(all 16) and every labelling, the point O = image of ∞ (which must lie on the circle of every line of L, including
O ∈ A when S is a 5-line) does not exist: Gröbner basis (1) in Q(ζ_10)[w, u, v, ρ] after saturating
ρ(ρ² − 1) Π(u − z_p), for finite O (u, v = O, conj O independent) and for O = ∞ (`runs/analyse_two5_all.log`).
In particular the theorem-only survivor (concentric pentagons with the 10 three-point lines
[[0,1,5],[0,4,8],[0,6,7],[1,3,7],[1,8,9],[2,3,9],[2,4,6],[2,7,8],[3,5,6],[4,5,9]], count 32) is not realisable.

### 4.3 Skeleton [5,5,5,5] (own 19)
13 (F, L) classes, count 32 each; every one has a derived Fano plane (hereditary SG violation) at some point,
and also orchard and Miquel violations.

## 5. Cited-table route (not needed)
With o(9) = 6, o(10) = 5 (Crowe–McKee 1968) and t3(10) = 12 (Burr–Grünbaum–Sloane 1974) the filters alone kill
36 of the 37 classes and the same single survivor remains (`runs/post_all.log`, "+cited" column); it is killed by
§4.2.  The cited values are therefore neither necessary nor sufficient.

## 6. Validation and cross-checks
* Skeleton list: identical to the independent list of `../verify_independent/n10/skeletons_check.py` (22).
* Lex symmetry breaking: unit test (all 0/1 vectors of length 5 with x <=_lex reverse(x)); skeleton [5] gives
  INFEASIBLE with and without it; skeleton [6] re-run without it (`runs/enum_1_nolex.log`, see addendum).
* Explicit t3(9) cap not load-bearing: skeletons [5], [6] give identical results without it (`runs/enum_*_nocap.*`);
  [5,5] disjoint: see addendum.
* Stage-2 line-set counts agree with an independent CP-SAT all-solution enumeration for every F-class with line
  sets (`runs/linesets_check.log`: 28, 28; 31, 30, 26; 32, 20).
* Agreement with the independent orderly enumeration `../n10-enum/orderly3` (SG caps only, no (8_3)/hereditary
  constraints in the search): their 9 candidates are the SG-only superset of mine — candidates 7, 8 = my F-classes
  0, 6 of [6] (28 line sets each), candidates 4, 5, 6 = my F-classes 19, 31, 58 of [5,5] disjoint (31, 30, 26 line
  sets), candidates 2, 3 = my skeleton 19; their candidates 0, 1 (skeletons [5,5]-2-common and [5,5,5]) contain an
  (8_3) and are excluded by my stage 1.  Also consistent with the v2 CP-SAT statuses (skeletons infeasible in v2
  have no (F, L) class here).
* `common.py` self-tests: Fano is SG-closed, (8_3) is not; a 7-point line is not a violation; antipodal 10-point
  structure passes all filters (count 33); cube structure is Miquel-closed, cube minus a block is not.
* Model validation on n = 8: not repeated here (done in v2 and verify_independent: optimum D + l = 39).

## 7. Files
`common.py` (utilities, filters), `skeletons_own.py` (+ `skeletons_own.json`, `skeletons_own.log`),
`lemmasABC.py`, `enum4.py` (two-stage enumerator), `postfilter.py`, `analyse_six_own.py`,
`involution_certificates.py`, `analyse_two5_own.py`, `two5_families.py`, `miquel_identity.py`, `t3_9_check.py`,
`linesets_cpsat_check.py`, `skel6_check.py`/`skel6_inspect.py` (the counter-example to the v2 hand proof),
`summarise.py`; outputs in `runs/` (`enum_<i>.json/.log`, `post_all.json/.log`, `post_7.*`,
`analyse_six_1.log`, `involution_certificates.log`, `analyse_two5_s4.log`, `analyse_two5_all.log`,
`two5_families.log`, `linesets_check.log`, `summary_table.md`).
CPU: stage 1 of [5,5] disjoint 367 s (1 worker); everything else seconds to a few minutes; robustness re-runs
(no symmetry breaking) up to ~1 h; machine shared (load 6–12).

## Addendum — robustness results (completed after the main text)
* Skeleton [6] without symmetry breaking: 9840 labelled solutions (1096 s), the same 8 F-classes and 8 (F, L)
  classes as the lex-leader run (30 labelled solutions, 4 s) — isomorphism-level comparison, 0 new classes.
* Skeleton [5,5] disjoint without the explicit t3(9) cap ((8_3) 8-subset constraints only): the same 59 F-classes
  and 16 (F, L) classes (`runs/enum_7_nocap.json`).  Together with [5] and [6] this shows the explicit cap is
  redundant in stage 1 (it is only used in the hand proof of §1.2, where it is justified by (STS9)+(8_3)).
* Numerical sanity check of the exact verdict for the theorem-only survivor (`pentagon_numeric.py`,
  `runs/pentagon_numeric.log`): random-start least squares over (ρ, O) finds zero-residual solutions only with
  ρ = ±1 (all ten points on one circle, |O| = 1), which is exactly the degenerate case removed by the saturation
  ρ(ρ² − 1) ≠ 0; with ρ restricted to [1.1, 50] or [0.02, 0.9] the best residual is >= 0.06 for every labelling.
  For O = ∞ the best residual is >= 0.30.
* No-lex re-runs of skeletons 5, 13, 19, 12 (`runs/enum_<i>_nolex.log`): 936 / 176 / 120 / 32 labelled solutions
  giving 6 / 10 / 9 / 2 F-classes and 0 / 0 / 13 / 0 (F, L) classes — identical to the lex-leader runs (30 / 34 / 18 / 8
  labelled solutions).  So every skeleton with solutions ([6], [5,5]-2-common, [5,5,5] x2, [5,5,5,5]) and the
  critical infeasible one ([5]) have been re-run without symmetry breaking with identical outcomes.
* Total CPU used by this agent: about 1 CPU-hour (main pipeline ~10 min; robustness re-runs ~40 min), on a machine
  shared with other agents (load 6–12).
