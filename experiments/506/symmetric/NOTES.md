# Erdős #506 — symmetric-orbit construction search (2026-09-05)

Strategy: configurations that are unions of orbits of a dihedral group D_m (m = 2..12):
regular m-gons on either mirror type (`A`: vertices at angles 2πk/m, `B`: at (2k+1)π/m),
generic 2m-point orbits (`G`, radius r and angle θ), the centre (`C`) and — in the Möbius
picture — the point at infinity (`I`).  1–3 orbits, continuous parameters (radii, angles),
n = 5..25.  Special parameter values are found as roots of coincidence conditions
(concyclic quadruples, collinear triples, quadruples containing the centre), one
representative per D_m-orbit of quadruples; the inversion centre O is optimised over all
pairwise intersections of blocks (circles(inv_O(P)) = |B(P)| − deg(O)).

## Files
- `symcore.py`   — point generation, exact-tolerance block computation, vectorised
                   inversion-centre search (`evaluate_fast`), inversion.
- `search.py`    — family enumeration, orbit representatives, vectorised 1-D root finding
                   along each free parameter from several "nice" base values, coordinate
                   iteration from roots (depth 2).  Output `results_m2_6.json`, `results_m7_12.json`
                   (+ logs `run_*.log`).
- `search2d.py`  — 2-parameter search: trace {condition A = 0} and find where a second
                   condition vanishes on it (`results_2d.json`).
- `verify.py`    — exact verification: parameters recognised exactly (nsimplify, or as an
                   algebraic root of the vanishing condition), coordinates in a sympy algebraic
                   number field (plus a quadratic extension for the inversion centre), exact
                   inversion, exact circle count.
- `independent_check.py` — second, independent exact check of the coordinate lists in
                   `records.json` (pure sympy radicals + simplify).
- `make_records.py` — merges results, verifies, writes `records.json`.
- `box17_check.py` — stand-alone exact check of the n = 8 record (rectangular box).

## Records (all verified exactly)

| n | circles | formula | configuration |
|---|---------|---------|---------------|
| 6 | **8**   | 9  | equilateral triangle + the 3 side midpoints (3 three-point lines, 3 four-point circles, 5 ordinary circles) |
| 7 | **11**  | 13 | equilateral triangle + 3 side midpoints + centroid (6 three-point lines = orchard maximum t3(7), 6 four-point circles, 5 ordinary circles) |
| 8 | **17**  | 19 | rectangular box a×b×c with c² = 3b² inscribed in a sphere, stereographically projected from the sphere point Q = (a, 2b, 0), which lies on the face plane x = a and on the two 3-vertex planes −x/a + y/b ± z/c = 1 (deg Q = 3, so 20 − 3 = 17). In the plane: two concentric similar rectangles (D_2-symmetric). Contains a 4-point line, two 3-point lines, 10 four-point circles, 6 ordinary circles. Segre's cube is the special case a = b = c, which only reaches 18 (max degree of a sphere point is 2 there). |

(Segre / the two-concentric-squares configuration gives 18 for n = 8; this is now beaten.)

Two-concentric-squares form of the n = 8 record (`concur.py` test): squares of radii 1 and ρ,
same orientation, reach 17 exactly at ρ ∈ {4 + √15, (4 + √7)/3} and their reciprocals
(box a = b, c² = 3a² resp. a² = 3c²; deg-3 point off the mirror axes); for every other ratio the
best is 18 (Segre's cube, a = b = c, is the ratio ρ = 2 + √3 and gives only 18).

## How the n = 8 record was found
The Möbius block structure of the cube (12 four-point planes + 8 three-point planes) is shared
by every rectangular box inscribed in a sphere.  For the box the sphere points of degree ≥ 3
were computed by hand: on the face x = a the two 3-vertex planes −x/a + y/b ± z/c = 1 meet the
sphere in a common point iff c² = 3b², giving Q = (a, 2b, 0).  No box has a sphere point of
degree 4 (checked: the only degree-3 points on a face are (a, ±2b, 0) [c² = 3b²] or
(a, 0, ±2c) [b² = 3c²], which cannot coexist).  The float search (`D2 G G`) rediscovers the
same configuration at many parameter values (e.g. r = 4, θ = π/6, and r = 4 + √15, θ = π/4).

## Methodological notes / pitfalls found on the way
- **Orbit-label bug (fixed).** The first version labelled points only by (type, index in orbit),
  so families with two orbits of the same type (A A, B B, G G, ...) got wrong "symmetry
  permutations"; the orbit reduction of conditions was ineffective (e.g. K = 7050 instead of
  698 for D7 A B B) and some conditions could be dropped.  Results produced by the buggy
  version are kept in `old_buggy/` for reference only; every family was rerun after the fix.
- **Quadruple coincidences are not the only events.**  The n = 8 record is a *concurrency*
  event: along the box family the block structure never changes, but at c² = 3b² three blocks
  pass through one point of the sphere.  Quadruple/triple conditions cannot see this; the 1-D
  search only found it because a "nice" base value (r = 4, θ = π/6) happens to be such a box.
  `concur.py` adds detection of these events (sign changes of the signed distance from the
  intersection points of block pairs to third blocks along a parameter slice).
- The inversion-centre candidate filter of the older `explore_families.py` accepted points of P;
  `symcore.best_centre` rejects candidates within 1e-5 of P, refines O from the actual points
  (double precision) so that the inverted set can be re-evaluated, and reports the blocks
  through O so that O can be recomputed exactly by `verify.py`.

## Relation to the combinatorial enumeration (`../theory/enum_n*_t*.json`)
- n = 8: the unique block structure with count ≤ 17 found by the theory enumeration is the
  12-block system AG(3,2) minus a parallel class with 3 lines {4-line, 3-line, 3-line} — exactly
  the structure realised by our box construction (`records.json`, n = 8).  So 17 is the minimum
  for that structure, and c(8) = 17 if the enumeration's necessary conditions are complete.
- n = 7: the unique structure with count ≤ 10 is the biplane 2-(7,4,2) (complements of the
  Fano lines; every point on 4 four-blocks).  It is **not Möbius-realisable**: invert about any
  point p; the 4 blocks through p become the 4 lines of a complete quadrilateral on the other 6
  points, and the 3 blocks missing p are the complements of the 3 Fano lines through p, i.e.
  precisely the three *diagonal quadruples* {opposite vertex pairs} of that quadrilateral — the
  Angle Lemma (`../theory/angle_lemma_check.py`) shows those cannot all be concyclic.  Our
  n = 7 record (11) realises 6 of the 7 four-blocks (all but {3 vertices, centroid}) — the best
  possible for this structure — so c(7) = 11 if the enumeration is complete.
  (`realise_design.py` was written to attack such realisation questions numerically; the run
  for the biplane produced no non-degenerate solution before its time limit, consistent with
  the proof above.)
- n = 6: the count-7 candidate needs the 3 diagonal quadruples of a complete quadrilateral to
  be concyclic — the same impossibility; our 8 is therefore optimal if the enumeration is
  complete.

## Negative results
(filled in below when the runs finish)
