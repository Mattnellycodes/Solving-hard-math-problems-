# Erdős #506 — results table, certificates and dependencies

Write-up of the session of 2026-09-05 (directory `experiments/506/`).  Every entry below is labelled with its
status and with the independent verification it rests on.  Nothing is called **proved** unless every
load-bearing step was reproduced by at least one implementation written from scratch by a different agent
(and, for n ≤ 9, adjudicated in `skeptics/ADJUDICATION.md`).  A companion research note is in
`docs/NOTE_506.md`.

**Convention** (Elliott's, also the Lean statement `erdos_506.variants.small_n` in
google-deepmind/formal-conjectures): P is a set of n ≥ 4 points of the plane, not contained in one line and
not contained in one circle.  A circle is *determined* by P if it passes through at least three points of P;
three collinear points determine no circle; lines with ≥ 3 points are allowed.  c(n) is the minimum number of
determined circles.  f(n) = C(n−1,2) + 1 − ⌊(n−1)/2⌋ is the Elliott/Purdy–Smith value, proved to equal c(n)
for n > 393.  Under a "no three collinear" convention the small values would differ (the extremal sets for
n = 6, 7, 8 contain collinear triples); nothing below is claimed under that convention.

## 1. The table

| n | c(n) | f(n) | status | extremal configuration (certificate in §2) | lower bound rests on (§3) | independent verification |
|---|---|---|---|---|---|---|
| 4 | **3** | 3 | proved | three collinear points + one point (or an antipodal pair + centre) | elementary (at most one collinear triple; two triples on one circle would make all four concyclic) | trivial; reproduced by every n = 4 enumeration |
| 5 | **5** | 5 | proved | four concyclic points + the centre, two diameters | hand proof (`theory/REPORT.md` §3.1) | adjudicator's enumeration `skeptics/adjudicator/enum67.py 5 5` (nothing ≤ 4) |
| 6 | **8** | 9 | proved (computer-assisted) — below f | triangle + feet of its altitudes | relaxation enumeration (SG-only caps) + Angle Lemma | `skeptics/ADJUDICATION.md` §2: CONFIRMED (two skeptics + adjudicator, own code) |
| 7 | **11** | 13 | proved (computer-assisted) — below f | orthocentric system + its three diagonal points | relaxation enumeration + Angle Lemma applied at a point | ADJUDICATION §3: CONFIRMED |
| 8 | **17** | 19 | proved (computer-assisted) — below f, and below Segre's projected cube (18) | inversion of the orthocentric configuration in a point where an altitude, the circumcircle and the nine-point circle concur | relaxation with Sylvester–Gallai-only caps (optimum D + ℓ = 39); no further geometry | ADJUDICATION §4: CONFIRMED (two CP-SAT formulations, a solver-free DFS, an orderly generation) |
| 9 | **25** | 25 | proved (computer-assisted) | 8 concyclic points in antipodal pairs + centre | relaxation (SG-only caps) + non-realisability of (8_3) + hereditary Sylvester–Gallai | ADJUDICATION §5: CONFIRMED; `verify_independent/cpsat_n9_sym.py` (OPTIMAL, 197 s) |
| 10 | **33** | 33 | proved (computer-assisted) | 9 concyclic points (4 antipodal pairs + 1) + centre | skeleton enumeration (SG-only caps) + (8_3) + hereditary SG + Miquel's theorem + exact Gröbner analysis of one surviving structure | `audit-n10-continue` (CONFIRMED, own SAT/DFS/Gröbner code); `audit-n10-cpsat` (independent pipeline, same conclusion); combinatorial part also by `n10-enum` + `audit-n10-enum` |
| 11 | ≤ 41 | 41 | best known; conjectured exact | antipodal (5 pairs + centre) | no non-trivial lower bound established | upper bound recounted exactly (three independent counters) |
| 12 | ≤ 51 | 51 | best known; conjectured exact | antipodal | — | idem |
| 13 | ≤ 61 | 61 | best known; conjectured exact | antipodal | — | idem |
| 14 | ≤ 73 | 73 | best known; conjectured exact | antipodal | — | idem |
| 15 | ≤ 85 | 85 | best known; conjectured exact | antipodal | — | idem |
| 16 | ≤ 99 | 99 | best known; conjectured exact | antipodal | — | idem |
| 17 … 393 | ≤ f(n) | f(n) | open | antipodal | — | — |
| ≥ 394 | f(n) | f(n) | proved (Elliott 1967; correction Purdy–Smith 2010) | antipodal | literature | — |

Sequence of proved values (n = 4, …, 10): **3, 5, 8, 11, 17, 25, 33**.  The exceptions to the formula
f(n) among 4 ≤ n ≤ 10 are exactly n = 6, 7, 8.  Conjecture (supported by the searches of §5 and by the
counting argument of `theory/REPORT.md` §5): c(n) = f(n) for all n ≥ 9.

Intermediate rigorous bounds at n = 10 that use fewer ingredients (purely combinatorial, no algebra):
c(10) ≥ 31 from Sylvester–Gallai + (8_3) alone, and c(10) ≥ 32 with Kelly–Moser's theorem o(m) ≥ 3m/7
(`n10-enum/orderly3`, CONFIRMED by `audit-n10-enum/v2-cpsat`).  Exactly two abstract structures with 32
circles survive every combinatorial test; both are the concentric-pentagon family, so geometry is
necessary for the last step.

## 2. Certificates (upper bounds)

All coordinate lists were recounted on 2026-09-05 by the write-up agent with `circles_exact.py`
(exact rational arithmetic; script `recount.py` in the session scratchpad, output reproduced here), in
addition to the counters of the agents that found them (`theory/certify.py`, `search-inversion/exact_verify.py`,
`search-local/verify_exact.py`) and of their auditors (`audit-search-inversion/v2/count506.py`,
`audit-search-local/count_exact.py`).  "circle sizes" = numbers of points on the determined circles;
"lines" = lines with ≥ 3 points.  None of the sets is contained in one line or one circle.

| n | circles | points | structure |
|---|---|---|---|
| 4 | 3 | (0,0), (1,0), (2,0), (0,1) | 3 ordinary circles; one 3-point line |
| 5 | 5 | (0,0), (1,0), (−1,0), (0,1), (0,−1) | one 4-point circle, 4 ordinary circles; two 3-point lines |
| 6 | 8 | (0,0), (20,0), (5,15), (5,0), (2,6), (10,10) | 3 four-point circles + 5 ordinary; 3 three-point lines (triangle (0,0),(20,0),(5,15) with the feet of its altitudes (5,0),(2,6),(10,10)) |
| 7 | 11 | (0,0), (20,0), (5,15), (5,0), (2,6), (10,10), (5,5) | 6 four-point circles + 5 ordinary; 6 three-point lines (the 7th point is the orthocentre) |
| 8 | 17 | (0,0), (0,5), (0,10), (0,15), (3,6), (3,9), (5,5), (5,10) | 11 four-point circles + 6 ordinary; one 4-point line (x = 0) and two 3-point lines {(0,0),(3,6),(5,10)}, {(0,15),(3,9),(5,5)} |
| 8 | 18 (Segre) | (±1,0), (0,±1), (±2,0), (0,±2) | 10 four-point circles + 8 ordinary; two 4-point lines — two concentric squares, Möbius-equivalent to a square prism on a sphere (the cube for ratio 2+√3) |
| 9 | 25 | (0,0), (25,0), (−25,0), (0,25), (0,−25), (7,24), (−7,−24), (24,7), (−24,−7) | one 8-point circle + 24 ordinary; 4 three-point lines (diameters) |
| 9 | 25 | (0,1), (0,3), (1,0), (1,4), (2,2), (3,0), (3,4), (4,1), (4,3) | same structure inside the 5 × 5 grid (unique 9-subset of the grid with 25 circles, `search-local`, confirmed by `audit-search-local`) |
| 10 | 33 | the nine points above + (15,20) | one 9-point circle + 32 ordinary; 4 three-point lines |
| 11 | 41 | + (−15,−20) | one 10-point circle + 40 ordinary; 5 lines |
| 12 | 51 | + (20,15) | one 11-point circle + 50 ordinary; 5 lines |
| 13 | 61 | + (−20,−15) | one 12-point circle + 60 ordinary; 6 lines |
| 14 | 73 | + (−7,24) | one 13-point circle + 72 ordinary; 6 lines |
| 15 | 85 | + (7,−24) | one 14-point circle + 84 ordinary; 7 lines |
| 16 | 99 | + (−24,7) | one 15-point circle + 98 ordinary; 7 lines |

Origin of the n = 8 record: triangle A = (0,3), B = (1,0), C = (3,0) with orthocentre H = (0,−1) and feet
D = (6/5,−3/5), E = (2,1), F = (0,0); the Möbius set {A,B,C,H,D,E,F,∞} has the *cube structure* (12
four-point blocks); the point O = (0,1) lies on three of its blocks (the altitude line AHF∞, the circumcircle
ABC and the nine-point circle DEF), and inverting in O gives the 8 points above (scaled by 10).  The whole
1-parameter family A = (0,3), B = (b,0), C = (3/b,0), 0 < b < √3, gives 17 for every b tried
(`theory/family17.py`: b = 1, 1/2, 3/2, 2/3, 5/4, 1/3, 7/5; not proved for all b).

General fact behind n ≥ 9: n − 1 points on a circle, ⌊(n−1)/2⌋ of the pairs collinear with a point x off
the circle (antipodal pairs and the centre), give 1 + C(n−1,2) − ⌊(n−1)/2⌋ = f(n) circles for every
n ≥ 4 (Lemma A below shows this is optimal among sets with an (n−1)-point circle).

## 3. Lower bounds: what each depends on

Notation (Möbius reformulation, `theory/REPORT.md` §1): a *block* is a circle or line through ≥ 3 points of P;
every triple lies in exactly one block; two blocks share ≤ 2 points; two lines share ≤ 1 point;
circles(P) = C(n,3) − D − ℓ with D = Σ_{|B| ≥ 4} (C(|B|,3) − 1) over all blocks and ℓ = number of lines with
≥ 3 points.  Inverting about p ∈ P maps the blocks through p to the ≥ 3-point lines of a real (n−1)-point set
(the *derived structure* at p), so Sylvester–Gallai applies to it.

### 3.1 Ingredients (all theorems; each re-proved or re-verified in the session)

| tag | statement | where re-verified |
|---|---|---|
| (ID) | the identity circles = C(n,3) − D − ℓ and the intersection axioms | four skeptic counters + adjudicator `recount.py`; every construction and hundreds of random rational sets |
| (SG) | Sylvester–Gallai, o(m) ≥ 1: the ≥ 4-blocks through p cover ≤ C(n−1,2) − 1 pairs of P∖{p}; the lines of P cover ≤ C(n,2) − 1 pairs; hereditary form for every subset not inside one block through p | adjudicator §6 (ii), (viii); `n10-continue/common.py`, `audit-n10-continue/lib.py` self-tests |
| (AL) | Angle Lemma: the three diagonal quadruples of a complete quadrilateral are never all concyclic | synthetic directed-angle proof + three independent algebraic certificates (`skeptics/adjudicator/angle_lemma.py`, `skeptics/c6-audit`, `skeptics/lemmas-audit`) |
| (8_3) | every 8 triples on 8 points pairwise sharing ≤ 1 point form the Möbius–Kantor configuration, which has no realisation by 8 distinct real points having exactly these 3-point lines (closing equation c² − c + 1 = 0) | `verify_independent/mk_unrealisable2.py`, `theory/eight_three_light.py`, `n10-cpsat/v2/facts.py` F3, `audit-n10-cpsat/mk_check.py`, `audit-n10-enum/v2-cpsat/mk_check.py`, `audit-n10-continue/mk_real_check.py`, adjudicator `configs.py` |
| (STS9) | every 11 triples on 9 points pairwise sharing ≤ 1 point contain an (8_3); with (8_3): at most 10 three-point lines on 9 real points | `n10-continue/t3_9_check.py`, `n10-cpsat/v2/facts.py` F4, `audit-n10-continue/sts9_check.py` (10080 labelled packings) |
| (MIQ) | Miquel's theorem: 8 distinct points of the inversive plane labelled by a cube never have exactly five concyclic faces (cross-ratio identity: product of the six face cross-ratios = 1) | `n10-continue/miquel_identity.py`, `n10-cpsat/v2/facts.py` F5, `audit-n10-cpsat/miquel_check.py`, `audit-n10-continue/miquel_identity.py` |
| (A),(B),(C) | largest-block lemmas: a block of size n−1 gives ≥ f(n) circles (equality iff antipodal type, n ≥ 5); size n−2 gives ≥ m² − m + 1 − 3⌊m/2⌋; size m in general gives ≥ 1 + r·C(m,2) − m·r(r−1)/4 − ℓ (r = n − m) | adjudicator §6 (iv)–(vi), `skeptics/lemmas-audit`, `n10-continue/lemmasABC.py`; for n = 10 not even needed (`n10-enum` enumerated blocks of size 7, 8, 9 directly: nothing) |
| (PAR) | parity: 14 three-point lines on 10 points are impossible (elementary packing bound D(10,3,2) = 13) | hand proof `n10-continue/REPORT.md` §1.2, checked by `audit-n10-continue` |

Cited values that are **not** load-bearing anywhere (used only in robustness re-runs or alternative kills):
Kelly–Moser o(m) ≥ 3m/7 (1958), Csima–Sawyer o(m) ≥ 6m/13 (1993), the table o(3..14) = 3,3,4,3,3,4,6,5,6,6,6,7
(Crowe–McKee 1968), and the orchard numbers t₃(3..12) = 1,1,2,4,6,7,10,12,16,19 (Burr–Grünbaum–Sloane 1974,
checked against OEIS A003035 by `audit-search-inversion/v2`).  Every enumeration below was run with
Sylvester–Gallai-only caps (o(m) = 1).

### 3.2 Per-n dependency

| n | argument | ingredients | programs (original / independent) |
|---|---|---|---|
| 5 | 4-block + point: ≤ 2 chords through the point; no 4-block: ≤ 2 lines (partial STS(5)) | (ID) | `theory/REPORT.md` §3.1 / `skeptics/adjudicator/enum67.py 5 5` |
| 6 | ≤ 7 circles forces three 4-blocks = complements of a perfect matching + 4 transversal 3-lines (a complete quadrilateral whose diagonal quadruples are concyclic) | (ID), (SG), (AL) | `theory/mobius_enum.py`, `verify_independent/cpsat_relaxation.py` / skeptics c6-*, adjudicator `enum67.py 6 7`, `cpsat_relax.py 6 sg` (optimum D + ℓ = 13) |
| 7 | ≤ 10 circles forces the (7,4,2) biplane of 4-blocks (ℓ = 7 would be a Fano plane of lines); at every point the derived structure is a complete quadrilateral with its three diagonal quadruples concyclic | (ID), (SG), (AL) | idem / skeptics c7-*, adjudicator `enum67.py 7 10`, `enum_big 7 m 10 1`, `cpsat_relax.py 7 sg` (27) |
| 8 | the relaxation optimum is D + ℓ = 39, i.e. 17 circles, already with SG-only caps; the unique count-17 structure is the cube structure with ∞ on 3 blocks, which is realised | (ID), (SG) | `theory/mobius_enum.py 8 17`, `theory/cpsat_model.py 8`, `verify_independent/cpsat_relaxation.py` / skeptics c8-*, adjudicator `cpsat_relax.py 8 sg` (OPTIMAL 39, 84 s), `enum_big 8 m 17 1` |
| 9 | ≤ 24 circles needs D + ℓ ≥ 60; blocks of size 6, 7, 8 excluded (Lemmas or direct enumeration); no 5-block ⇒ a point in ≥ 8 four-blocks ⇒ derived (8_3); with a 5-block the unique structure (two 5-blocks through a point + 12 four-blocks, count 24) has a derived Fano plane at every degree-7 point, contradicting hereditary SG | (ID), (SG), (8_3) | `theory/mobius_enum2.py 9 24`, `theory/hereditary_sg.py`, `verify_independent/cpsat_n9_sym.py` / skeptics c9-*, adjudicator `enum_big 9 m 24 {0,1}` (m = 4…8, with and without caps), `postproc.py`, `check_m4.py` |
| 10 | see §3.3 | (ID), (SG), (8_3), (STS9), (MIQ), (PAR), exact algebra; (A)–(C) optional | `n10-continue/` / `audit-n10-continue/`, `audit-n10-cpsat/`; combinatorial part `n10-enum/orderly3/` / `audit-n10-enum/`, `audit-n10-enum/v2-cpsat/` |

### 3.3 n = 10 in detail (c(10) ≥ 33)

Reduction: ≤ 32 circles ⇔ D + ℓ ≥ 88; (SG) gives ℓ ≤ 14 and per-point caps (35 pairs on the derived 9-set,
44 on the inverted 10-set Q_p = (P∖p) ∪ {∞'}).

1. **Largest block ≥ 7** — excluded: size 9 gives ≥ 33 (Lemma A), size 8 ≥ 45 (Lemma B; ≥ 39 by the cruder
   count with ℓ ≤ 14), size 7 ≥ 40 (Lemma C).
   Independently, `n10-enum` enumerated the families with a block of size 7, 8 or 9 under SG-only caps:
   none reaches D + ℓ ≥ 88.
2. **All blocks of size ≤ 4** — (STS9) gives ≤ 10 four-blocks through a point, so b₄ = 25, every point in
   exactly 10 four-blocks, ℓ ≥ 13; (SG) on Q_p bounds the lines through each point, and (PAR) finishes
   (hand proof `n10-continue/REPORT.md` §1.2, checked by reasoning in `audit-n10-continue`).  Independently:
   `n10-enum/orderly3/all4.py` enumerates the case (1 class, 27 four-blocks, killed by (8_3));
   `audit-n10-enum/v2-cpsat/all4.py` (through a maximum-degree point and all PSTS(9) classes: same single class);
   `audit-n10-cpsat` (4 families with b₄ = 25, none admits ≥ 13 lines).
3. **Skeletons** (families of blocks of size ≥ 5, pairwise sharing ≤ 2 points, under the SG cap): 22
   isomorphism classes, computed by four independent programs (`verify_independent/n10/skeletons_check.py`,
   `n10-continue/skeletons_own.py`, `audit-n10-continue/skeletons.py`, `audit-n10-enum/v2-cpsat/pipe.py`);
   18 have largest block 5 or 6.  For each: stage 1 = complete enumeration of the 4-block families
   (constraints: triples in ≤ 1 block, SG cap, (8_3) on every 8-subset of every derived 9-set, D ≥ 88 − ℓ_max);
   stage 2 = all line sets with D + ℓ ≥ 88; reduction to isomorphism classes.  Results, identical in
   `n10-continue` (CP-SAT + own DFS), `audit-n10-continue` (pysat CaDiCaL/Glucose + own DFS, full-orbit
   blocking) and `audit-n10-cpsat` (CP-SAT with own lex-leader symmetry breaking):

   | skeleton (block sizes; pairwise intersections) | 4-block families (classes) | (F,L) classes with D+ℓ ≥ 88 | killed by |
   |---|---|---|---|
   | [5] | infeasible | 0 | — |
   | [6] | 8 | 8 (counts 31/32, ℓ = 12/13) | Miquel closure (all 8); also exact involution analysis of the 6-circle (`n10-continue/analyse_six_own.py`: both families with line sets are unrealisable) |
   | [5,5] sharing 2 points | 6 | 0 | no line set reaches ℓ ≥ 88 − D |
   | [5,5] sharing 1 point | infeasible | 0 | — |
   | [5,5] disjoint | 59 (6 with b₄ = 20, 53 with b₄ = 19) | 16 (counts 31/32) | Miquel closure (14), (8_3)+orchard cap (1), **1 survivor → exact algebra** |
   | [5,6] (two skeletons) | infeasible | 0 | — |
   | [6,6] | 1 | 0 | no line set |
   | [5,5,5] (three skeletons) | 2, 10, infeasible | 0 | no line set |
   | [5,5,6] | infeasible | 0 | — |
   | [5,5,5,5] (|Aut| = 128) | 9 | 13 (count 32) | hereditary SG (a derived Fano plane at some point); also orchard/Miquel |
   | [5,5,5,5] (two skeletons, |Aut| = 8) | 1, infeasible | 0 | no line set |
   | [5,5,5,6], [5⁵], [5⁶] | infeasible | 0 | — |

   Total: 37 (F,L) classes; 36 fail a theorem-level necessary condition; 1 survives all of them.
4. **The survivor** (the same isomorphism class in all three pipelines and in `n10-enum`'s list —
   checked again by the write-up agent with VF2; |Aut(F)| = 40, |Aut(F,L)| = 20):
   two disjoint 5-blocks S = {0,1,2,3,4}, R = {5,6,7,8,9}, the 20 four-blocks
   `[0,1,5,6],[0,1,7,8],[0,2,5,7],[0,2,6,9],[0,3,5,8],[0,3,7,9],[0,4,5,9],[0,4,6,8],[1,2,6,7],[1,2,8,9],[1,3,5,9],[1,3,6,8],[1,4,5,7],[1,4,6,9],[2,3,5,6],[2,3,7,8],[2,4,5,8],[2,4,7,9],[3,4,6,7],[3,4,8,9]`
   (every point in 9 rich blocks, D = 78) and the 10 three-point lines
   `[0,1,9],[0,2,8],[0,6,7],[1,3,7],[1,5,8],[2,4,6],[2,5,9],[3,4,5],[3,6,9],[4,7,8]` (ℓ = 10, count 32;
   the only other line set in its Aut(F)-orbit is `[1,2,5],[0,3,6],[4,5,6],[0,4,7],[3,5,7],[1,4,8],[2,6,8],[2,3,9],[1,7,9],[0,8,9]`).
   *Exact exclusion.*  After a Möbius transformation S and R lie on distinct circles; the three normal forms
   (concentric circles, parallel lines, two lines through a point) turn the 20 concyclicity conditions
   {s,s',r,r'} into an integer linear system whose rational kernel is the constants and whose torsion is
   Z/10, so the only realisations of the 4-block family are two concentric regular pentagons (angles
   2πk/10, aligned or rotated by π/5; 8 labellings), at an arbitrary radius ratio ρ ≠ 0, ±1 — confirmed
   numerically (20 four-point circles, isomorphic block structure, `verify_independent/n10/pentagon_check2.py`
   and the write-up agent's check).  The 10 lines must be circles through the image O of ∞; the resulting
   polynomial system in (O, ρ), saturated by ρ(ρ² − 1) ≠ 0, has Gröbner basis {1} for finite O and for
   O = ∞, for all 8 labellings, in three independent formulations: `n10-continue/analyse_two5_own.py`
   (Q(ζ₁₀)[w,u,v,ρ]), `audit-n10-cpsat/two5_exact.py` (Q[u,v,ρ,z]/(Φ₂₀(z)); in fact all 16 (F,L) classes of
   the skeleton are so excluded) and `audit-n10-continue/survivor_exact_real.py` (real coordinates with
   c = cos π/5, s = sin π/5 as ring variables; positive controls return non-unit ideals).  Numerical
   corroboration: `verify_independent/n10/pentagon_scan.py`, `n10-continue/pentagon_numeric.py`,
   `audit-n10-continue/survivor_numeric.py` (zero residual only at the degenerate ρ = ±1).
5. **Soundness remark.** A real configuration whose incidences refine one of the enumerated structures is
   itself one of the enumerated structures (the enumeration is over all structures with count ≤ 32), and
   every exclusion used is a necessary condition, so the case analysis is exhaustive.

Robustness checks recorded: stage 1 re-run without symmetry breaking for [5], [6], [5,5]-2-common, [5,5,5] ×2,
[5,5,5,5] (identical classes; e.g. [6]: 9840 labelled solutions = sum of orbit sizes); without the explicit
t₃(9) cap (identical); DFS line-set counts = CP-SAT counts; the SG-only candidate list of `n10-enum` (9 classes)
is the superset of the (8_3)-filtered lists; the antipodal 33-circle structure passes every filter (positive
control); the (n = 8, 9) analogues of every pipeline reproduce the known results.

### 3.4 Defects found in agent reports (do not cite these passages)

* `n10-cpsat/v2/REPORT.md` §4.0: the hand proof for the single-5-block skeleton is wrong (a slack of 2 is
  overlooked; a 22-four-block family with the stated degree caps exists).  The conclusion is true by
  computation (stage 1 infeasible in all three pipelines, with and without symmetry breaking).
* The argument "all-4-block case: ℓ ≤ 12 by t₃(10) = 12" (task context, `theory/REPORT.md` §6) overlooks 12
  three-point lines plus one 4-point line; replaced by item 2 above.
* `n10-continue/REPORT.md` §4.2 prints the survivor with a mis-transcribed block list (four printed lines lie
  inside printed blocks); the machine record `n10-continue/runs/post_all.json` #12 is correct and is the
  structure printed above (isomorphic to the records of both audits).
* `theory/angle_lemma_check.py` omits the chart with a line perpendicular to L1; complete certificates are in
  `skeptics/adjudicator/angle_lemma.py`, `skeptics/c6-audit`, `skeptics/lemmas-audit`.
* `theory/cpsat_hereditary.py` encodes hereditary SG without the hypothesis "subset not inside one block
  through p" and is unsound (declares the antipodal 9-point structure infeasible); never used for a result.
* `theory/REPORT.md` §3.5: Lemma C for n = 9, m = 6 gives 27 (26 with SG-only caps), not 28.
* `explore_families.py` (float exploration) accepted inversion centres coinciding with points of P until fixed;
  no result depends on it.

## 4. Structure theorems (by-products)

| statement | status |
|---|---|
| Every 8-point set with ≤ 18 circles carries the cube structure (12 four-point blocks, SQS(8) minus a parallel class); the sets with 17 circles are exactly the cube-structure sets with ∞ on three blocks | computed (`theory/enum_n8_t18.json`); supported by the adjudicator's independent enumeration (2 classes with count ≤ 18, the second needing an (8_3) at infinity) and by `audit-n10-enum/v2-cpsat` (n = 8, target 18: same two classes) |
| Every 9-point set with 25 circles consists of 8 concyclic points and a 9th point lying on 4 lines through pairs of them (the chords need not be diameters) | computed (`theory/enum2_n9_t25.json`); independently reproduced by `n10-enum/orderly3` (n = 9, target 25, all block sizes: 9 candidates, only the 8-block structure survives (8_3)/hereditary SG) and by `audit-n10-enum/v2-cpsat` (same classes and line-set counts) |
| Every 10-point set with ≤ 32 circles: none (§3.3); no classification of the 33-circle sets was attempted | — |

## 5. Construction searches (n = 9 … 16): negative results

Two independent searches (`search-inversion`, `search-local`; audits `audit-search-inversion`,
`audit-search-local`, both CONFIRMED) looked for n-point sets with fewer than f(n) circles by subset selection
plus inversion-centre optimisation over coincidence-rich universes (integer lattice spheres x² + y² + z² = N,
grids and grids with inversions, triangular lattices, regular polygons at several radii, orthocentric
closures, Platonic/Archimedean solids, conics, polygon-diagonal closures, coincidence closures; 68 universes,
580 simulated-annealing runs; one-parameter family scans detecting concurrency/concyclicity events, which
rediscover the n = 8 record).  Outcome: nothing below f(n) for any 9 ≤ n ≤ 16; every set attaining f(n) has
the antipodal Möbius type (one (n−1)-block + ⌊(n−1)/2⌋ lines); the best genuinely different structures are far
above f(n) (n = 10: 37 from two concentric pentagons at the pentagram ratio; 39 on the N = 5 lattice sphere).
Exhaustive results: all 2,042,975 nine-subsets of the 5 × 5 grid with the inversion centre ranging over all
arrangement points have minimum 25 (85 subsets, all antipodal type, unique Euclidean minimiser) and nothing in
26..31; all 3,268,760 ten-subsets have minimum 45 (`audit-search-local/grid5_n10.py`).  The lattice sphere
N = 2 (cuboctahedron) gives 31, 44, 58 for n = 9, 10, 11 (exhaustive, reproduced).  These are search results,
not proofs.

## 6. Reproduction

Upper bounds: `python3 experiments/506/circles_exact.py` (built-in examples) or `analyze(points)` on any list
above.  Lower bounds n ≤ 9: the command list in `skeptics/ADJUDICATION.md` §8.  n = 10: `n10-continue/`
(`skeletons_own.py`, `enum4.py`, `postfilter.py`, `analyse_two5_own.py`; ≈ 1 CPU-hour) and the audits
`audit-n10-continue/` (`stage1_sat.py`, `stage2_lines.py`, `postfilter.py`, `survivor_exact_real.py`; the last
takes ≈ 15 min) and `audit-n10-cpsat/`.  Combinatorial bound c(10) ≥ 32: `n10-enum/orderly3/oe.py`,
`all4.py`, `hcheck.py` (≈ 17 CPU-minutes).  Environment: Python 3.11 with ortools 9.15, python-sat, sympy 1.14,
numpy, networkx (`docs/ENVIRONMENT.md`).
