# Adjudication — Erdős #506, claims of `experiments/506/theory/REPORT.md`

Adjudicator, 2026-09-05.  Inputs: the theory report and its code, the ten skeptic verdicts
(two lenses × five claims, directories `skeptics/c6-*`, `c7-*`, `c8-*`, `c9-*`, `lemmas-*`), and my
own from-scratch verification code in `skeptics/adjudicator/` (nothing imported from the theory
agent or from any skeptic).  Standard applied: a claim is **CONFIRMED** only if every load-bearing
step is supported by an independent computation or proof (a skeptic's *and*, for the load-bearing
steps, my own); otherwise **GAP** (with what is missing) or **REFUTED**.

Convention adjudicated (the task's convention, = Elliott's): c(n) = minimum number of circles
containing ≥ 3 points of an n-point planar set that is not contained in one line or one circle;
collinear triples determine no circle; lines with ≥ 3 points are allowed.  Under a "no three
collinear" convention the small values would differ; that is outside the claims.

## 0. Verdict table

| claim | statement | verdict |
|---|---|---|
| c6 | c(6) = 8 (construction 8 circles; no 6-set with ≤ 7) | **CONFIRMED** |
| c7 | c(7) = 11 | **CONFIRMED** |
| c8 | c(8) = 17 (relaxation optimum D+ℓ = 39; unique count-17 structure = cube with 3 lines) | **CONFIRMED** |
| c9 | c(9) = 25 | **CONFIRMED** |
| lemmas | (i) identity, (ii) derived/hereditary SG caps, (iii) Angle Lemma, (iv) Lemma A, (v) Lemma B, (vi) Lemma C, (vii) (8_3) uniqueness + non-realisability, (viii) hereditary SG | **CONFIRMED** (with recorded, non-load-bearing defects in two of the theory's *certificate programs* and two typos; see §6) |

Overall: the values c(5)=5, c(6)=8, c(7)=11, c(8)=17, c(9)=25 stand.  No counterexample was found
by any of the ten skeptics' construction searches (≈ 2·10⁹ exact grid/rational subsets for n=6,7,8,9
plus structured/numerical searches) and none by me.  The only outstanding items are cosmetic fixes
to the theory write-up and its (non-load-bearing) certificate code.

## 1. What the whole argument rests on (checked line by line)

1. **Möbius identity** circles(P) = C(n,3) − D − ℓ, D = Σ_{|B|≥4}(C(|B|,3)−1) over all blocks
   (circles and lines), ℓ = #lines with ≥ 3 points.  Every triple lies in exactly one block, so
   #blocks = C(n,3) − D; circles = #blocks − ℓ.  Verified exactly on every construction and on
   hundreds of random rational sets by four skeptics and by my `recount.py`.
2. **Necessary conditions** (the relaxation): (C1) two blocks share ≤ 2 points; (C3) two lines
   share ≤ 1 point; a line of size ≥ 4 is a block; a 3-point line is not inside a larger block;
   (C2) Sylvester–Gallai caps with o = 1 only: at each p the ≥4-blocks through p cover ≤ C(n−1,2)−1
   pairs (inversion about p turns them into the ≥3-point lines of a non-collinear real set, which
   has an ordinary line); lines cover ≤ C(n,2)−1 pairs.  All are theorems; the relaxation is a
   superset of the realisable structures, so every real configuration with count ≤ target has its
   exact structure among the enumerated survivors.  The Kelly–Moser/Csima–Sawyer values quoted
   from memory in the report are **not load-bearing**: every enumeration below was (re)run with
   o = 1 only and gives the same survivors.
3. **Exhaustive enumeration** of the relaxation at each n (three to four mutually independent
   programs per n, see §2–§5), leaving a short list of survivors.
4. **Non-realisability of each survivor** by a theorem: Angle Lemma (n=6,7), Fano plane at a point
   or at infinity ⇒ no ordinary line ⇒ contradicts Sylvester–Gallai (n=7 lines, n=9), (8_3) at a
   point (n=9, no-5-block case).  Each of these was re-proved independently (§5).

## 2. Claim c6 — CONFIRMED

Skeptics (c6-construct, c6-audit): own exact recounts (8 circles, 3 lines), own exhaustive
enumerations (Python brute force over all 82 block families and all line sets, and a z3 model),
own Angle-Lemma certificates (slope algebra A1−A2 = 2(m4−m3)(1+m1m2); a two-chart Gröbner
certificate closing the theory's chart gap), exact brute force over 1.0·10⁹ six-subsets of a 7×7
grid and a 97-point rational set (minimum 8, never 7).

My own checks (`adjudicator/`):
* `recount.py` (blocks identified purely by point sets from integer 3×3/4×4 determinants):
  n=6 set → `circles=8 lines=3 D=9 identity=True degenerate=False`.
* `enum67.py 6 7` (plain DFS over all families of 4-/5-subsets pairwise sharing ≤ 2 points, exact
  ℓ_max, canonical forms by all 720 permutations): regimes `none` and `sg` both give
  **exactly one** isomorphism class with count ≤ 7: blocks {0123},{0145},{2345}, ℓ_max = 4, count 7
  (15 labelled) — i.e. a complete quadrilateral with all three diagonal quadruples on blocks.  A
  5-block never reaches 7, so Lemma A is not even needed.  Control `enum67.py 5 5`: only the
  4-block + 2 lines structure (count 5), nothing ≤ 4, so c(5) = 5.
* `cpsat_relax.py 6 sg`: OPTIMAL max D+ℓ = 13 (min 7).
* `angle_lemma.py`: **single complete chart** L1: y = 0, L_i: x = m_i y + c_i (covers every line
  not parallel to L1, vertical included; c2 = 0, c3 = 1 by translation/scaling).  The three
  concyclicity determinants factor as (non-degeneracy factors) × E_k with
  E1 = m2m3 − m2m4 − m3m4 − 1, E2 = m2m3 − m2m4 + m3m4 + 1, E3 = m2m3 + m2m4 − m3m4 + 1, and
  `Groebner(Q1,Q2,Q3, t*nondeg-1) = [1]` — the Angle Lemma holds for every complete
  quadrilateral (this closes the chart gap of the theory's own `angle_lemma_check.py`).  Any two
  of the conditions force a perpendicularity (lex bases contain m2 = 0 and m3m4+1, etc.), exactly
  the synthetic proof's mechanism; the two-condition family is realisable (Q3 = 78 ≠ 0), so the
  lemma is sharp.
* The reduction: in the unique candidate each point is on exactly two of the four 3-lines and any
  two lines share exactly one point, so the four lines are pairwise non-parallel and no three are
  concurrent; the three 4-blocks cannot be lines (a 4-point line would share 2 points with a
  transversal 3-line, violating (C3)), hence they are circles and are exactly the three diagonal
  quadruples (my `configs.py` (E)-type check and both skeptics).  Angle Lemma ⇒ impossible ⇒
  c(6) ≥ 8; the construction gives 8.

## 3. Claim c7 — CONFIRMED

Skeptics (c7-construct, c7-audit): own recounts (11), own labelled enumeration of all 6149 block
families and an own CP-SAT (max D+ℓ = 28 no caps / 27 with SG caps), own Gröbner/factorisation
proofs of the Angle Lemma in independent parametrisations (including the vertical-line charts and
the collinear-instead-of-concyclic degenerations), Fano non-realisability (det = 2), exhaustive
grid searches (5×5, 6×6, 7×7, 8×6, 9×5, 10×4: nothing ≤ 11 on grids, minimum 13).

My own checks:
* `recount.py`: n=7 set → `circles=11 lines=6 D=18 identity=True`.
* `enum67.py 7 10`: regime `none`: 2 classes with count ≤ 10 — the (7,4,2) biplane (7 four-blocks,
  ℓ_max = 7, count 7) and "biplane minus a block" (6 four-blocks, ℓ_max = 7, count 10); regime
  `sg`: **exactly one** class — the biplane with ℓ_max = 6, count 8.  ℓ = 7 on 7 points means seven
  pairwise ≤1-sharing lines covering all 21 pairs, i.e. a Fano plane of straight lines: no
  ordinary line, contradicting Sylvester–Gallai (`configs.py` (C): unique 7-triple system, det 2).
* `enum_big 7 4 10 1` / `7 5 10 1` / `7 6 10 1` (own C enumerator, largest block forced): 6 / 0 / 0
  labelled survivors, the 6 being the biplane copies containing {0,1,2,3} (7!/168·7/35 = 6 ✓).
* `cpsat_relax.py 7 sg` → OPTIMAL 27 (min 8); `7 none` → 28 (min 7).
* `configs.py` (E): in the biplane, at **every** point p the four blocks through p give four
  3-point derived lines pairwise sharing exactly one point with every other point on exactly two
  of them (complete quadrilateral), and the three blocks avoiding p are exactly the three diagonal
  quadruples.  Non-parallel/non-concurrent as in §2; blocks avoiding p (circles or lines) invert
  to genuine circles, blocks through p invert to lines, so the argument is independent of which
  blocks are lines.  Angle Lemma (my certificate) ⇒ impossible ⇒ c(7) ≥ 11.

## 4. Claim c8 — CONFIRMED

Skeptics (c8-construct, c8-audit): own recounts (17: 11 four-point + 6 three-point circles, lines
of sizes 4,3,3), two own CP-SAT formulations (all OPTIMAL: max D+ℓ = 39 ⇒ 17 under SG-only and
under strong caps; 44 ⇒ 12 with no caps), an own solver-free DFS (393,328 nodes, 36 labelled
count-17 families = one class, the cube) and an own orderly generation validated on n = 5,6,7
(cube unique at ≤ 17; cube + one 10-block/8-line structure at ≤ 18), exact grid searches
(6×6 … 8×6, 7×7, 11×4: minimum 18 both Euclidean and after optimal inversion), polyhedral and
augmentation searches (best 17), least-squares realisers (control finds the cube).

My own checks:
* `recount.py`: n=8 set → `circles=17 lines=3 D=36 identity=True`, blocks of size ≥ 4: 12 (all
  degrees 6), lines {0123} (4 points), {0,4,7}, {3,5,6}.  Control: two concentric squares → 18.
* `cpsat_relax.py 8 sg` (triple-free formulation, x/y variables, pairwise conflicts, SG caps 20/27):
  `status=OPTIMAL max D+ell=39 bound=39 => min circles=17`, 84 s, 2 workers; `8 none`: OPTIMAL 44
  (min 12, SQS(8) + 2 lines, all degrees 7) — so the SG cap is load-bearing and valid.
* `enum_big 8 m 17 1` for m = 4,5,6,7: found 36 / 0 / 0 / 0; the 36 are one class (VF2), the cube
  (D = 36, ℓ_max = 3, count 17, 36 = 210·12/70 labelled copies containing {0,1,2,3} ✓).  Since
  target 17 includes every count ≤ 16, **no relaxation structure has ≤ 16 circles**.
  `enum_big 8 4 18 1`: 216 = 36 (cube) + 180 (10 four-blocks, all degrees 5, D = 30, ℓ_max = 8,
  count 18), identical to the skeptics' two classes.
* The cube structure with ∞ on three blocks is realised by the construction, so 17 is attained.

## 5. Claim c9 — CONFIRMED

Skeptics (c9-construct, c9-audit): own recounts (25, two witnesses), two independent C enumerators
with the largest block fixed, run under strong, SG-only and **no** caps (m = 8,7,6: nothing;
m = 5: exactly 30 labelled = one class; m = 4: 6480 (caps) / 55440 (no caps) labelled families with
count ≤ 24, every one with a point in 8 four-blocks), own hereditary-SG checkers (derived Fano at
all eight degree-7 points), own (8_3) proofs (frame + Gröbner, closure argument, 81-chart analysis),
own verification of Lemmas A/B/C (hand + exact random sets), exhaustive 5×5-grid 9-subsets
(minimum 25, attained only by the antipodal configuration), random grid/closure/annealing searches
(best 25–26).

My own checks:
* `recount.py`: antipodal set → 25; a second witness (8 unit-circle points in pairs collinear with
  (1/3, 0), i.e. four non-diameter chords) → 25.
* `enum_big 9 m 24 c` (own C enumerator; validated above on n = 7, 8 where it reproduces the
  known classes and labelled counts exactly):
  - m = 8, 7, 6 with SG caps and with no caps: `found=0` (m = 6: 120,970 nodes, 0 families with
    D ≥ Dmin).  This settles blocks of size 6, 7, 8 without Lemmas A, B, C.
  - m = 5: SG caps: 34,760,422 nodes, 390 families checked, **30 found**; no caps: 41,748,880
    nodes, 72,360 checked, 30 found.  `postproc.py` (VF2): **one class** — two 5-blocks sharing a
    point + twelve 4-blocks, degrees (7⁸, 2), D = 54, own ℓ_max = 6, count 24 — the report's
    candidate.  My `configs.py` (D) and `postproc.py`: at each of the 8 degree-7 points p, with q
    the degree-2 point and S = P∖{p,q}, the derived lines restricted to S are 7 three-point lines
    covering all 21 pairs (a Fano plane) and S lies in no block through p; inverting about p gives
    7 real non-collinear points with no ordinary line — contradiction with Sylvester–Gallai.
  - m = 4 (no 5-block), SG caps: 714,105,490 nodes, `found=6480`, histogram
    (D, ℓ, count, maxdeg) = {(51,11,22,8): 6120, (54,11,19,8): 360} — identical to c9-audit's
    SG-only run; **every** family has a point in exactly 8 four-blocks (`check_m4.log`).  This is
    also forced by counting: count ≤ 24 ⇒ D + ℓ ≥ 60, ℓ ≤ ⌊36/3⌋ = 12 ⇒ b₄ ≥ 16 ⇒ Σdeg ≥ 64 > 63.
    At such a point the derived structure is 8 triples on 8 points pairwise sharing ≤ 1 point
    (9 impossible: `configs.py` (A) finds 0), all degrees 3; with only 4-blocks through p its
    collinear triples are exactly these 8, so the frame {0,1,2,5} is in general position and the
    remaining incidences force c² − c + 1 = 0 (`configs.py` (B): Gröbner [a+c−1, b+c−1, c²−c+1],
    no real root; frame-free least squares never finds a non-collinear real solution).
  - Consistency: `enum_big 9 8 25 1` returns the antipodal structure (D = 55, ℓ = 4, count 25),
    `9 5 25 1` returns only the count-24 candidate.
* Hence every relaxation structure with ≤ 24 circles is non-realisable and c(9) = 25.

## 6. Claim "lemmas" — CONFIRMED (with recorded defects)

| lemma | status | independent support |
|---|---|---|
| (i) identity | proved; exact | 4 skeptic counters + mine on all constructions and random sets |
| (ii) derived SG caps, hereditary form | proved (inversion argument re-derived by 3 skeptics and by me); SG-only caps suffice everywhere | all enumerations rerun with o = 1; n = 8 and n = 9 also with **no** caps where relevant |
| (iii) Angle Lemma | proved: synthetic directed-angle proof checked by 4 skeptics; 3 independent algebraic certificates (c6-audit two charts, c7-construct/c7-audit/lemmas-audit parametrisations, my single complete chart) | unit ideal in every chart; sharp |
| (iv) Lemma A | proved (hand, re-derived here); brute-force confirmed n = 5..8 | not load-bearing for n ≤ 9 (enumerations cover all block sizes); "equality iff" clause false only at n = 4 |
| (v) Lemma B | proved (I re-derived ℓ ≤ m by the three cases: no xy-line; xy-line with 1 point of B; 4-point line {x,y,a,b}); brute-force 0 violations (n = 6,7,8) | not load-bearing for n ≤ 9 |
| (vi) Lemma C | proved (chunk argument re-derived here); 0 violations over 17,362 abstract structures and 399+ exact random sets, tight cases exist | not load-bearing for n ≤ 9; §3.5 quotes 28, correct value 1+36−10 = 27 (still > 24) |
| (vii) (8_3) | 840 labelled 8-triple systems = one class = Möbius–Kantor (3 skeptics + me); non-real: c² − c + 1 = 0 (3 independent frames/charts + 81-chart analysis + mine); every realisation by distinct non-collinear points has exactly the 8 lines (mk_propagation closure argument) | load-bearing for c(9), m = 4 case, and for the n = 8 "≤ 18 ⇒ cube" by-product |
| (viii) hereditary SG | immediate from (ii) applied to the inverted subset; hypothesis "S not inside one block through p" is essential | load-bearing for c(9), m = 5 case; the theory's `hereditary_sg.py` includes the hypothesis (sound) |

Defects found (none affects the five value claims):
1. `theory/angle_lemma_check.py` uses L1: y = 0 and finite slopes for L2..L4, so it cannot
   represent a line perpendicular to L1 — precisely the case through which the synthetic proof
   passes.  Its "unit ideal" is a chart artefact for the two-condition statement (in the
   perpendicular chart the two-condition family is a real 2-parameter family).  The lemma itself
   is true; complete certificates exist (c6-audit chart A+B, lemmas-audit chart A+B, my
   `adjudicator/angle_lemma.py`).  Fix: replace the certificate by one of these.
2. `theory/cpsat_hereditary.py` encodes (H1)/(H2) without the hypothesis "S is not contained in
   one block through p" (resp. one line).  I read the code: `sum(terms) <= 6` is imposed for every
   7-subset S and every Fano labelling, so an 8-block containing p ∪ S violates it.  The
   lemmas-audit demonstration (`hereditary_bug_demo.py`) shows the realisable antipodal 9-point
   structure is declared INFEASIBLE.  The tool was not used for n ≤ 9 but is proposed for n = 10,
   where it would wrongly exclude every structure with a block of size ≥ 9.  Fix: add
   7·x_{S∪{p}} (resp. 8·x_{S∪{p}} for (H2)) to the right-hand side, or skip subsets inside a block.
3. REPORT §3.5: Lemma C bound for n = 9, m = 6 is 27 (26 with SG-only caps), not 28.
4. REPORT Lemma 2.5: the equality clause "iff B is a circle …" fails for n = 4 (three collinear
   points + one point also gives f(4) = 3); true for n ≥ 5.
5. Lemma 2.4's sentence "8 real points span at most 7 three-point lines" is justified in the
   report only for exact structures; that is all its uses need (the general statement is the
   classical orchard value t₃(8) = 7).

## 7. Residual caveats (not gaps in the claims)

* Convention: the claims hold under the stated convention (points not all on one line or circle;
  collinear triples determine no circle); the erdosproblems.com wording is ambiguous, as noted by
  c6-audit.
* "New record" / novelty statements were not adjudicated (literature access was limited).
* The by-product structure theorems (every 8-set with ≤ 18 circles carries the cube structure;
  every 25-circle 9-set is 8 concyclic points + a point on 4 chords) were not part of the
  adjudicated claims.  The n = 8 one is supported by my enumeration (2 classes at ≤ 18, the second
  needing an (8_3) at infinity); the n = 9 one would need the target-25 enumeration for largest
  block 4, 6, 7 as well, which I did not run.
* Several skeptic runs were still executing when they reported (numerical realisers, extra grid
  pools, a z3 cross-check, one sympy Gröbner run); none of them is load-bearing, and each such
  step is covered by a finished independent computation listed above.

## 8. Reproduction (my code; all in `experiments/506/skeptics/adjudicator/`)

```
python3 recount.py                       # exact recounts: 8, 11, 17, 18(control), 25, 25
python3 angle_lemma.py                   # single-chart certificate, unit ideal (9 min under load)
python3 configs.py                       # (8_3) uniqueness/non-realisability, Fano, derived Fano at n=9, biplane pattern
python3 enum67.py 5 5 ; enum67.py 6 7 ; enum67.py 7 10
gcc -O2 -o enum_big enum_big.c
./enum_big 7 4 10 1 ; ./enum_big 8 4 17 1 ; ./enum_big 8 4 18 1        # validation
./enum_big 9 8 24 1 ; ./enum_big 9 7 24 1 ; ./enum_big 9 6 24 1 ; ./enum_big 9 5 24 1 ; ./enum_big 9 5 24 0 ; ./enum_big 9 4 24 1
python3 postproc.py n9_m5_t24_c1.txt ; python3 check_m4.py n9_m4_t24_c1.txt
python3 cpsat_relax.py 6 sg ; 7 sg ; 7 none ; 8 sg 900 ; 8 none 600   # 13, 27, 28, 39, 44
```
Logs: `recount.log`, `angle_lemma.log`, `configs.log`, `enum67.log`, `val_*.txt/.err`,
`n9_*.txt/.err`, `postproc_*.log`, `check_m4.log`, `cpsat_chain.log`.
