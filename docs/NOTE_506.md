# The minimum number of circles determined by n points: exact values for n ≤ 10 (Erdős problem #506)

Draft research note, 2026-09-05.  Prepared in the repository `Solving-hard-math-problems-`
(directory `experiments/506/`).  All computations were carried out by AI agents (Claude) under a workflow in
which every load-bearing step was re-implemented from scratch by at least one independent agent and, for
n ≤ 9, adjudicated by a separate skeptic workflow; the write-up was checked against the logs.  **A human
review is still required before this is posted or submitted.**  Bibliographic details of the older papers
were taken from citing papers and search snippets (the originals of Elliott 1967 and Bálintová–Bálint 1994
were not accessible from the session) and should be checked against the originals.

## Abstract

Let c(n) be the minimum number of circles passing through at least three points of an n-point set in the
plane that is contained neither in a line nor in a circle (three collinear points determine no circle).
Elliott (1967), with the correction of Purdy and Smith (2010), proved c(n) = f(n) := C(n−1,2) + 1 − ⌊(n−1)/2⌋
for n > 393; no exact value was recorded for 4 ≤ n ≤ 393.  We determine

    c(4), …, c(10) = 3, 5, 8, 11, 17, 25, 33.

The values c(6) = 8, c(7) = 11 and c(8) = 17 lie below f(n) = 9, 13, 19 (and c(8) = 17 is below the 18 circles of
Segre's projected cube); c(9) = f(9) and c(10) = f(10).  The extremal configurations for n = 6, 7, 8 come from a
triangle with the feet of its altitudes and its orthocentre; for n = 8 the record set is an inversion of the
orthocentric configuration in a point where an altitude, the circumcircle and the nine-point circle concur.
The lower bounds are computer-assisted: the incidence structure of circles and lines of a point set is a
"Möbius block structure" (every triple in exactly one block, two blocks share at most two points), inversion
about a point turns the blocks through it into the lines of a real point set, and Sylvester–Gallai-type
constraints plus a finite list of non-realisable sub-configurations (the Fano plane, the Möbius–Kantor
configuration (8_3), the "Angle Lemma" for complete quadrilaterals, Miquel's theorem) cut the enumeration
down to a handful of abstract structures, each of which is then excluded by an exact argument.  For n = 10
a single abstract structure with 32 circles survives all necessary conditions — two concentric regular
pentagons with ten three-point lines — and is excluded by a Gröbner-basis computation.  Explicit integer
coordinates certify every upper bound; extensive searches found nothing below f(n) for 11 ≤ n ≤ 16.  We
conjecture c(n) = f(n) for all n ≥ 9.

## 1. The problem

**Problem (Erdős 1961; erdosproblems.com #506).**  Let P be a set of n points in the plane, not all on one
line and not all on one circle.  A circle is *determined* by P if it passes through at least three points
of P.  What is the minimum number c(n) of circles determined by such a set?

Conventions.  Collinear triples determine no circle, and lines through three or more points of P are
allowed; this is the convention of Elliott's theorem and of the Lean statement
`erdos_506.variants.small_n` in google-deepmind/formal-conjectures (the erdosproblems.com text is ambiguous
between this and "no three collinear").  Under "no three collinear" the small values would be different —
the extremal sets below for n = 6, 7, 8 contain collinear triples — and nothing in this note concerns that
variant.  The sequence starts at n = 4 (three points are always concyclic or collinear).

Known results.  Elliott [El67] proved that for n > 393 the points determine at least C(n−1,2) circles;
Purdy and Smith [PS10] observed that the correct bound is

    f(n) = C(n−1,2) + 1 − ⌊(n−1)/2⌋,

attained by n − 1 points on a circle in antipodal pairs together with the centre (the ⌊(n−1)/2⌋ triples
{centre, x, −x} are collinear and determine no circle), and that Elliott's proof gives it for n ≥ 394.
Segre observed (reported by Elliott) that the projection of a cube gives an 8-point set with fewer than
C(7,2) = 21 circles.  Related asymptotic results concern ordinary circles (circles through exactly three
points): Elliott, Bálintová–Bálint [BB94], Zhang [Zh11], Lin–Makhul–Mojarrad–Naslund–Swanepoel [LMMNS17],
Lin–Makhul–Mojarrad–Schicho–Swanepoel–de Zeeuw [LMMSSZ18], all for n large.  To our knowledge no exact value
of c(n) for any 4 ≤ n ≤ 393 appears in the literature, and no OEIS entry exists (teorth/erdosproblems issue
#297 asks for exactly this computation).

## 2. Results

**Theorem 1.**  With the convention above,

| n | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| c(n) | 3 | 5 | **8** | **11** | **17** | 25 | 33 |
| f(n) | 3 | 5 | 9 | 13 | 19 | 25 | 33 |

In particular the formula f(n) fails exactly for n = 6, 7, 8 within 4 ≤ n ≤ 10, and the projected cube
(18 circles) is not optimal for n = 8.

**Theorem 2 (structure).**  (a) Every 8-point set determining at most 18 circles has the *cube structure*:
its 4-point circles-or-lines form the Steiner quadruple system SQS(8) minus a parallel class (12 blocks);
the sets with exactly 17 circles are the cube-structure sets in which three of the twelve blocks are
straight lines.  (b) Every 9-point set determining 25 circles consists of 8 concyclic points and a ninth
point lying on four lines each through two of the eight (the chords need not be diameters and the point need
not be the centre).  (c) Every 10-point set determining at most 32 circles: there is none (Theorem 1).

**Corollary.**  The sequence c(4), c(5), … begins 3, 5, 8, 11, 17, 25, 33.

**Conjecture 3.**  c(n) = f(n) for all n ≥ 9, i.e. the exceptions to the Elliott–Purdy–Smith formula are
exactly n = 6, 7, 8.

Status of the proofs.  Upper bounds are explicit integer configurations, recounted in exact arithmetic by
three or more independent programs.  Lower bounds are computer-assisted case analyses (§5–§6); every
enumeration and every algebraic certificate was reproduced by an independent implementation (§7).  For
n ≥ 11 only c(n) ≤ f(n) is known; the value f(n) was attained but never beaten in the searches of §8.

## 3. Constructions

All coordinates are integers; the counts were verified exactly (rational arithmetic) by
`experiments/506/circles_exact.py` and by the independent counters listed in §7.  "Ordinary" = circle through
exactly three points.

**3.1 (n = 4, 5).**  (0,0), (1,0), (2,0), (0,1): three circles.  (0,0), (±1,0), (0,±1): one 4-point circle and
four ordinary circles, five in all (two diameters through the centre).

**3.2 (n = 6, eight circles).**  A = (0,0), B = (20,0), C = (5,15) and the feet of the altitudes
D = (5,0), E = (2,6), F = (10,10).  The three circles with diameters AB, BC, CA each contain two feet
(three 4-point circles); the three sides are 3-point lines; the remaining triples give 5 ordinary circles
(20 = 3·4 + 3·1 + 5 triples), so the set determines 3 + 5 = 8 circles, while f(6) = 9.

**3.3 (n = 7, eleven circles).**  The six points of 3.2 and the orthocentre H = (5,5).  The complete
quadrangle A, B, C, H has six sides (3-point lines through the diagonal points D, E, F) and the six circles
with diameters XY (X, Y ∈ {A,B,C,H}) each contain two of D, E, F: 6 four-point circles + 5 ordinary
circles = 11, while f(7) = 13.

**3.4 (n = 8, seventeen circles).**  (0,0), (0,5), (0,10), (0,15), (3,6), (3,9), (5,5), (5,10):
eleven 4-point circles and six ordinary circles; the lines with ≥ 3 points are x = 0 (four points),
{(0,0),(3,6),(5,10)} and {(0,15),(3,9),(5,5)}.  Origin: take the triangle A = (0,3), B = (1,0), C = (3,0),
its orthocentre H = (0,−1) and the feet D = (6/5,−3/5), E = (2,1), F = (0,0).  In the inversive plane the
eight points A, B, C, H, D, E, F, ∞ carry the cube structure (twelve 4-point blocks: the six lines of the
quadrangle ABCH together with ∞, and the six diameter circles).  The point O = (0,1) is the reflection of H
in BC (so it lies on the circumcircle ABC), the midpoint of AH (so it lies on the nine-point circle DEF), and
lies on the altitude AHF; inverting in O turns these three blocks into lines and gives the set above (after
scaling by 10).  The same works for A = (0,3), B = (b,0), C = (3/b,0) for every 0 < b < √3 tried
(b = 1, 1/2, 3/2, 2/3, 5/4, 1/3, 7/5), so the record is attained by a one-parameter family.
For comparison, Segre's projected cube — two concentric squares (±1,0), (0,±1), (±2,0), (0,±2), or the cube
itself projected along a 4-fold axis, ratio 2 + √3 — gives 18 circles, and 18 is the minimum over the whole
Möbius class of the cube (every non-vertex point of the circumsphere lies on at most two of the cube's
twenty 3-point planes; exact computation, `experiments/506/literature/segre_cube_check.py`).

**3.5 (n = 9, 10: the antipodal sets).**  (0,0), (25,0), (−25,0), (0,25), (0,−25), (7,24), (−7,−24), (24,7),
(−24,−7): one 8-point circle, four 3-point lines (diameters), 24 ordinary circles, 25 in all.  Adding
(15,20) gives a 10-point set with 33 circles.  Another 25-circle 9-set lies in the 5 × 5 grid:
(0,1), (0,3), (1,0), (1,4), (2,2), (3,0), (3,4), (4,1), (4,3) (eight lattice points on (x−2)² + (y−2)² = 5 plus the
centre); it is the unique 9-subset of the grid with 25 circles (exhaustive search, §8).

**3.6 (n = 11 … 16).**  Continuing with (−15,−20), (20,15), (−20,−15), (−7,24), (7,−24), (−24,7) on
x² + y² = 625 gives 41, 51, 61, 73, 85, 99 = f(11), …, f(16) circles.  In general n − 1 points on a circle,
⌊(n−1)/2⌋ pairs of them collinear with a point x off the circle, give exactly f(n) circles for every n ≥ 4.

## 4. The Möbius reformulation and the lemmas

Work in the inversive plane (plane plus ∞).  For a finite set P not contained in one circle or line, a
**block** is a circle or line containing at least three points of P, regarded as a subset of P.

**4.1 Identity.**  Every triple of P lies in exactly one block; two blocks share at most two points; two lines
share at most one point.  Writing D = Σ_{|B| ≥ 4} (C(|B|,3) − 1) over all blocks and ℓ for the number of lines
(blocks through ∞),

    circles(P) = C(n,3) − D − ℓ.

So c(n) = min C(n,3) − D − ℓ over all *realisable* pairs (F, L), where F is the family of blocks of size ≥ 4
and L the family of lines.  Equivalently, c(n) = |B(P)| − ℓ where B(P) is the set of all blocks; inverting P
about a point O ∉ P replaces ℓ by the number of blocks through O.

**4.2 Derived structures and Sylvester–Gallai.**  Inverting about p ∈ P maps the blocks through p to the
lines with ≥ 3 points of the real (n−1)-point set P∖{p} (the *derived structure* at p) and the other blocks to
circles.  Since P∖{p} is not collinear, the Sylvester–Gallai theorem gives it an ordinary line, so the blocks of
size ≥ 4 through p cover at most C(n−1,2) − 1 pairs of P∖{p}; likewise the lines of P cover at most C(n,2) − 1
pairs of P, and the same holds for the inverted 10-set Q_p = (P∖{p}) ∪ {∞'} (image of P ∪ {∞}).  *Hereditary
form:* every subset S of the derived set that is not contained in one line has an ordinary line of its own.
Consequences: no derived structure contains a Fano plane (7 three-point lines on 7 points cover all pairs), and
a point of P lies in at most t₃(n−1) four-point blocks, where t₃(m) is the maximum number of 3-point lines of
m real points.  (Quantitative bounds — Kelly–Moser o(m) ≥ 3m/7, Csima–Sawyer o(m) ≥ 6m/13, the Crowe–McKee
table of o(m), the orchard numbers t₃(m) of Burr–Grünbaum–Sloane — are *not* used in any proof below;
o(m) ≥ 1 suffices throughout.)

**4.3 Angle Lemma.**  Let L₁, …, L₄ be four lines, pairwise non-parallel, no three concurrent, with vertices
P_ij = L_i ∩ L_j.  The three "diagonal quadruples" {P₁₂,P₃₄,P₁₃,P₂₄}, {P₁₂,P₃₄,P₁₄,P₂₃}, {P₁₃,P₂₄,P₁₄,P₂₃} are
not all concyclic.  *Proof.*  With directed angles mod π and θ_i the direction of L_i, four points A, B, C, D
are concyclic iff ∡(CA,CB) = ∡(DA,DB); the three conditions read θ₁ + θ₄ ≡ θ₂ + θ₃, θ₁ + θ₃ ≡ θ₂ + θ₄,
θ₁ + θ₂ ≡ θ₃ + θ₄ (mod π).  The first two give 2θ₁ ≡ 2θ₂, i.e. L₁ ⊥ L₂; the first and third give L₁ ⊥ L₃;
hence L₂ ∥ L₃, a contradiction.  ∎  (Algebraic certificates: in a complete chart the three concyclicity
determinants factor into non-degeneracy factors times E₁ = m₂m₃ − m₂m₄ − m₃m₄ − 1, E₂ = m₂m₃ − m₂m₄ + m₃m₄ + 1,
E₃ = m₂m₃ + m₂m₄ − m₃m₄ + 1, and the ideal saturated by the non-degeneracy conditions is the unit ideal;
`experiments/506/skeptics/adjudicator/angle_lemma.py`.)

**4.4 The Möbius–Kantor configuration.**  Every family of 8 triples on 8 points pairwise sharing at most one
point is isomorphic to the configuration (8_3) = {i, i+1, i+3 mod 8} (all 840 labelled families form one
S₈-orbit; |Aut| = 48), and (8_3) has no realisation by 8 distinct real points whose 3-point lines are exactly
its 8 lines: with a projective frame on four points containing no line of the configuration the remaining
incidences force the closing equation c² − c + 1 = 0.  (Classical — Möbius 1828, Kantor 1881; see Grünbaum
[Gr09] — and re-proved from scratch in the session, `experiments/506/verify_independent/mk_unrealisable2.py` and
four further independent implementations.)  Consequences: 8 real points span at most 7 three-point lines
(t₃(8) = 7); no point of a 9-point set lies in 8 four-point blocks; and, since every 11 triples on 9 points
pairwise sharing ≤ 1 point contain an (8_3) (exhaustive check of the 10080 labelled packings, three
independent programs), t₃(9) ≤ 10, so no point of a 10-point set lies in more than 10 four-point blocks.

**4.5 Largest-block lemmas.**  Let B be a largest block, |B| = m, r = n − m.
(A) If m = n − 1 then circles(P) ≥ f(n), with equality (n ≥ 5) iff B is a circle and the remaining point lies
on ⌊(n−1)/2⌋ lines through pairs of B.  (The C(n−1,2) triples {x,a,b} lie in distinct blocks; if B is a line no
other line exists; if B is a circle the lines through x cut B in disjoint pairs.)
(B) If m = n − 2 then circles(P) ≥ m² − m + 1 − 3⌊m/2⌋, which exceeds f(n) for n ≥ 7.
(C) (chunk lemma) The number of blocks meeting B in exactly two points is at least r·C(m,2) − m·r(r−1)/4;
hence circles(P) ≥ 1 + r·C(m,2) − m·r(r−1)/4 − ℓ.  (For a ≠ b in B the other blocks through {a,b} partition
P∖B into "chunks"; for fixed a the chunks of different b share ≤ 1 point, so the pairs inside chunks are
counted at most once.)
For n = 10 these give: a block of size 9 forces ≥ 33 circles, size 8 forces ≥ 45, size 7 forces ≥ 40.

**4.6 Miquel closure.**  If 8 distinct points of the inversive plane are labelled by the vertices of a cube
and five of the six faces are concyclic quadruples, so is the sixth (Miquel's six-circle theorem; the
cross-ratio identity behind it — the product of the six face cross-ratios equals 1 — was re-verified
symbolically).  Every exact block structure must therefore be Miquel-closed.

**4.7 Elementary packing facts.**  Fourteen 3-point lines on 10 points are impossible (they would cover 42 of
the 45 pairs, and the three uncovered pairs cannot give all ten points an odd "leave" degree); equivalently the
packing number D(10,3,2) = 13.

## 5. Lower bounds for n ≤ 9 (proof sketches)

The method: enumerate all abstract pairs (F, L) satisfying 4.1 and the Sylvester–Gallai caps of 4.2 with
C(n,3) − D − ℓ ≤ target (a *relaxation* — a superset of the realisable structures), then exclude each
survivor by a theorem.  Because a real configuration whose incidences refine a listed structure is itself a
listed structure, and every exclusion is a necessary condition, this is exhaustive.

*n = 5.*  A 4-block B and a point x: the six triples {x,a,b} lie in distinct blocks; if B is a line there is no
other line (6 circles), if B is a circle at most two chords through x are lines (≥ 5 circles).  No 4-block:
10 blocks and at most two lines (a partial triple system on 5 points), ≥ 8.  So c(5) = 5.

*n = 6.*  circles = 20 − 3q − ℓ with q ≤ 3 four-blocks (a 5-block gives ≥ 9 by (A)) and ℓ ≤ 4.  Seven circles
forces q = 3, ℓ = 4: the three 4-blocks are the complements of a perfect matching and the four lines are the
transversal triples, i.e. a complete quadrilateral whose three diagonal quadruples are concyclic —
impossible by the Angle Lemma.  So c(6) = 8.

*n = 7.*  Ten circles needs D + ℓ ≥ 25 with ℓ ≤ 6 (seven lines would form a Fano plane, contradicting
Sylvester–Gallai), so D ≥ 19 with no block of size ≥ 5 (by (A), (B)): seven 4-blocks whose complements are
the lines of a Fano plane.  At any point p the four blocks through p invert to a complete quadrilateral on the
other six points and the three blocks avoiding p are its diagonal quadruples — impossible by the Angle Lemma.
So c(7) = 11.

*n = 8.*  The relaxation alone suffices: the maximum of D + ℓ under 4.1 and the Sylvester–Gallai caps (with
o = 1) is 39, i.e. at least 17 circles (CP-SAT optimality proof, reproduced by an orderly enumeration and a
solver-free search).  The unique structure with 17 circles is the cube structure with three blocks through
∞, realised in 3.4.  Every structure with ≤ 18 circles is either the cube structure or a 10-block structure
requiring an (8_3) of lines, which proves Theorem 2(a).

*n = 9.*  Twenty-four circles needs D + ℓ ≥ 60 with ℓ ≤ 11.  Blocks of size 6, 7, 8 are excluded ((A)–(C), or
directly by enumeration).  Without a 5-block, D = 3b₄ ≥ 49 forces some point into 8 four-blocks, whose
derived structure is an (8_3): impossible (4.4).  With a 5-block the enumeration leaves exactly one
structure (two 5-blocks through a common point plus twelve 4-blocks; 24 circles), and at each of its eight
points of degree 7 the derived structure contains a Fano plane on the other seven points, contradicting the
hereditary Sylvester–Gallai property.  So c(9) = 25.  The same enumeration with target 25 (all block sizes,
including the all-4-block case) leaves only the 8-block-plus-4-lines structure, proving Theorem 2(b).

## 6. The lower bound c(10) ≥ 33 (proof sketch)

Thirty-two circles or fewer means D + ℓ ≥ 88; Sylvester–Gallai gives ℓ ≤ 14 and per-point caps (35 pairs on
each derived 9-set, 44 on each inverted 10-set Q_p and on the lines of P).

**Step 1 (largest block ≥ 7).**  Excluded by 4.5; independently, the enumeration of all block families with a
block of size 7, 8 or 9 under the caps finds nothing with D + ℓ ≥ 88.

**Step 2 (all blocks of size ≤ 4).**  By 4.4 each point lies in ≤ 10 four-blocks, so b₄ ≤ 25; with D = 3b₄
and ℓ ≤ 14 we need b₄ = 25, every point in exactly 10 four-blocks, and ℓ ≥ 13.  Sylvester–Gallai on Q_p
shows that each point lies on at most 4 lines of P, hence at most one 4-point line exists, and in either case
some point p has four 3-point lines through it, so Q_p carries 14 three-point lines — impossible by 4.7.
(Independently: direct enumeration of the case under the caps gives a single family of 27 four-blocks, killed
by (8_3); a third implementation finds four families with b₄ = 25, none admitting 13 lines.)

**Step 3 (skeletons).**  The families of blocks of size ≥ 5 on 10 points, pairwise sharing ≤ 2 points and
satisfying the caps, fall into 22 isomorphism classes (four independent programs), 18 of which have largest
block 5 or 6: [5]; [6]; [5,5] with 0, 1 or 2 common points; [5,6] ×2; [6,6]; [5,5,5] ×3; [5,5,6]; [5,5,5,5] ×3;
[5,5,5,6]; [5⁵]; [5⁶].  For each skeleton all completions by 4-blocks were enumerated (triples in ≤ 1 block,
the caps, the (8_3) condition on every 8-subset of every derived 9-set, D ≥ 74), then all line sets with
D + ℓ ≥ 88, and the pairs (F, L) were reduced to isomorphism classes.  Three independent pipelines (CP-SAT
with lex-leader symmetry breaking; a SAT solver with full-orbit blocking; CP-SAT with a different symmetry
breaking, plus an orderly generation for the combinatorial part) agree class by class:

| skeleton | 4-block families | (F,L) classes with ≤ 32 circles | excluded by |
|---|---|---|---|
| [5]; [5,5] sharing 1; [5,6] ×2; [5,5,5] (one of three); [5,5,6]; [5,5,5,5] (one of three); [5,5,5,6]; [5⁵]; [5⁶] | none | 0 | — |
| [5,5] sharing 2; [6,6]; [5,5,5] (two); [5,5,5,5] (one) | 6; 1; 2, 10; 1 | 0 | no line set is large enough |
| [6] | 8 | 8 | Miquel closure (each; also an exact analysis of the involutions cut on the 6-circle) |
| [5,5,5,5] with |Aut| = 128 | 9 | 13 | hereditary Sylvester–Gallai (a derived Fano plane) |
| [5,5] disjoint | 59 | 16 | Miquel closure (14), (8_3) + t₃(9) ≤ 10 (1), **one survivor** |

**Step 4 (the survivor).**  The single structure passing every necessary condition is, up to isomorphism,
S = {0,1,2,3,4} and R = {5,6,7,8,9} as 5-point circles, the twenty 4-point circles

    {0,1,5,6} {0,1,7,8} {0,2,5,7} {0,2,6,9} {0,3,5,8} {0,3,7,9} {0,4,5,9} {0,4,6,8} {1,2,6,7} {1,2,8,9}
    {1,3,5,9} {1,3,6,8} {1,4,5,7} {1,4,6,9} {2,3,5,6} {2,3,7,8} {2,4,5,8} {2,4,7,9} {3,4,6,7} {3,4,8,9}

(every point on nine of the 22 rich blocks, D = 78) and the ten 3-point lines

    {0,1,9} {0,2,8} {0,6,7} {1,3,7} {1,5,8} {2,4,6} {2,5,9} {3,4,5} {3,6,9} {4,7,8}

(ℓ = 10; 120 − 78 − 10 = 32 circles; |Aut(F)| = 40, |Aut(F,L)| = 20).  The 4-block family is realisable: after a
Möbius transformation S and R lie on distinct circles, and in each of the three normal forms (concentric
circles |z| = 1, |z| = ρ; parallel lines; two lines through a point) the condition "{s,s',r,r'} concyclic" is a
linear relation on angles (α_s + α_s' ≡ β_r + β_r' mod 2π) or on coordinates.  The twenty relations form an
integer system whose rational kernel is the constants and whose torsion group is Z/10 (Smith normal form),
so the only realisations are two concentric regular pentagons — vertex angles 2πk/10, aligned or rotated by
π/5, eight labellings — at an arbitrary radius ratio ρ ≠ 0, ±1.  The ten lines of L must then be circles
through the point O that ∞ is mapped to.  For each labelling the ten conditions "O on the circle through the
three points of the line" are polynomial in (O, ρ) over Q(ζ₁₀); after saturating by ρ(ρ² − 1) ≠ 0 the Gröbner
basis is {1}, both for finite O and for O = ∞.  Hence the structure is not realisable, and c(10) ≥ 33.  The
Gröbner computation was done in three independent formulations (complex variables over Q(ζ₁₀); over
Q[u,v,ρ,z]/(Φ₂₀(z)); real coordinates with cos(π/5), sin(π/5) as ring variables, with positive controls
returning non-unit ideals) and corroborated numerically (least squares finds zero residual only at ρ = ±1,
where all ten points are concyclic).  In fact all 16 (F,L) classes of the disjoint-5-block skeleton, and both
6-block families that admit line sets, are also excluded by exact algebra, so Miquel's theorem is not the only
route.

Purely combinatorially (no algebra): Sylvester–Gallai and (8_3) alone give c(10) ≥ 31; adding Kelly–Moser's
o(m) ≥ 3m/7 gives c(10) ≥ 32, and exactly two abstract structures with 32 circles survive every
combinatorial test — both are the concentric-pentagon family above with one of its two line-set orbits — so
geometry is needed, and suffices, for the last step.

## 7. Verification and dependencies

What the proofs depend on: the identity 4.1 and the intersection axioms; the Sylvester–Gallai theorem
(only o(m) ≥ 1, in derived and hereditary form); the Angle Lemma (n = 6, 7); the non-realisability of (8_3)
(n = 9, 10) and the derived fact t₃(9) ≤ 10 (n = 10); Miquel's theorem (n = 10); the elementary packing facts
4.7 and, optionally, the largest-block lemmas 4.5; exhaustive enumerations of finite relaxations (orderly
generation, CP-SAT, SAT); Gröbner-basis computations (sympy) for n = 10.  No cited numerical table
(Kelly–Moser, Csima–Sawyer, Crowe–McKee, Burr–Grünbaum–Sloane) is load-bearing; the enumerations were also
run with those values and give the same survivors.

Independence.  For each n = 6, 7, 8, 9 the enumeration was performed by at least three programs written
independently (orderly generation `theory/mobius_enum*.py`; CP-SAT `verify_independent/cpsat_relaxation.py`
and `cpsat_n9_sym.py`; the skeptics' and the adjudicator's C/Python/z3 enumerators in
`skeptics/`), and every lemma has at least two independent proofs or certificates (adjudication record
`skeptics/ADJUDICATION.md`).  For n = 10 the complete argument was produced by `n10-continue/` and reproduced
class by class by `audit-n10-continue/` (pysat/CaDiCaL instead of CP-SAT, own DFS, own Smith normal form, real
Gröbner bases) and by `audit-n10-cpsat/`; the combinatorial part was reproduced by the orderly generation
`n10-enum/orderly3/` and by `audit-n10-enum/` (own canonical labelling + CP-SAT).  The write-up agent
re-verified all coordinate certificates and the survivor's internal consistency and isomorphism class.
Several defects in intermediate agent reports were found and are listed in `experiments/506/RESULTS.md` §3.4;
none affects the results (in each case the flawed passage was replaced by a computation or a corrected
argument).

## 8. Searches for n = 11, …, 16

Two independent searches (`search-inversion/`, `search-local/`, both audited) looked for sets with fewer
than f(n) circles: subset selection with simulated annealing or exhaustive enumeration over coincidence-rich
"universes" (integer lattice spheres x² + y² + z² = N, grids and grids closed under inversions, triangular
lattices, regular polygons at several radii, orthocentric closures, Platonic and Archimedean solids, conics,
diagonal-intersection closures of polygons — 68 universes, 580 annealing runs), always optimising the
inversion centre over all pairwise intersections of blocks; and one-parameter family scans detecting
concurrency and concyclicity events (which rediscover the n = 8 record).  Result: nothing below f(n) for any
9 ≤ n ≤ 16; every set attaining f(n) has the antipodal Möbius type (one (n−1)-point block and ⌊(n−1)/2⌋
lines); the best structurally different sets are far above (n = 10: 37 circles from two concentric pentagons
at the pentagram ratio, 39 from the 24 points of x² + y² + z² = 5).  Exhaustively: over all 2,042,975 nine-subsets
of the 5 × 5 grid and all inversion centres the minimum is 25 (85 subsets, all antipodal type), with no value
in 26..31; over all 3,268,760 ten-subsets the minimum is 45.  A counting argument shows why exceptions to f(n)
should be confined to small n: a set beating f(n) needs D + ℓ ≥ C(n−1,3) + ⌊(n−1)/2⌋ with largest block at
most ≈ n/2, i.e. almost every derived structure must be a near-extremal orchard configuration, while blocks of
size 4 alone can supply at most about 3n·t₃(n−1)/4 + t₃(n) < C(n−1,3) for n ≥ 10.

## 9. Open problems

1. Prove Conjecture 3, or at least c(11) = 41 and c(12) = 51.  The n = 10 pipeline (skeletons of blocks of
   size ≥ 5, CP-SAT/SAT completion with (8_3) and hereditary Sylvester–Gallai constraints, Miquel closure,
   exact algebra for survivors) should extend to n = 11 (need D + ℓ ≥ 125; largest block ≤ 6) with a few
   CPU-hours; n = 12 needs stronger symmetry breaking or local-structure enumeration (derived structures
   must be near-extremal orchard configurations, which are classified for ≤ 12 points).
2. Classify the 10-point sets with exactly 33 circles (is the antipodal type the only one, as for n = 9?).
3. Lower the threshold 393 of Elliott's theorem by re-running his inversion argument with the modern
   ordinary-line bounds (Csima–Sawyer, Green–Tao) and the exact small values of c(n).
4. Is every 8-point set with 17 circles Möbius-equivalent to a member of the orthocentric family of 3.4?
   (Theorem 2(a) fixes the block structure; the realisation space of the cube structure with three
   concurrent blocks was not analysed.)
5. Determine c(n) under the "no three collinear" convention.

## 10. Draft OEIS entry

    %N  Minimum number of circles determined by n points in the plane that are not all on one line and not
        all on one circle; a circle is determined if it passes through at least three of the points
        (three collinear points determine no circle).
    %S  3, 5, 8, 11, 17, 25, 33
    %O  4,1
    %C  a(n) = binomial(n-1,2) + 1 - floor((n-1)/2) for n > 393 (Elliott 1967, with the correction of Purdy
        and Smith 2010); this value is attained for every n >= 4 by n-1 points on a circle in antipodal pairs
        together with the centre, so it is an upper bound for all n.  It equals a(n) for n = 4, 5, 9, 10 but
        not for n = 6, 7, 8, where a(n) = 8, 11, 17 < 9, 13, 19.  Conjecture: a(n) = binomial(n-1,2) + 1 -
        floor((n-1)/2) for all n >= 9.
    %C  The terms a(4)..a(10) were established by computer-assisted proofs (2026): the incidence structure
        of circles and lines is enumerated as a partial "Möbius block structure" constrained by the
        Sylvester–Gallai theorem applied after inversion about each point, and the surviving abstract
        structures are excluded by non-realisability arguments (Fano plane, Möbius–Kantor configuration,
        an angle lemma for complete quadrilaterals, Miquel's theorem, Gröbner bases).  Each step was
        reproduced by independent programs.
    %D  P. D. T. A. Elliott, On the number of circles determined by n points, Acta Math. Acad. Sci. Hungar.
        18 (1967), 181-188.
    %D  G. B. Purdy and J. W. Smith, Lines, circles, planes and spheres, Discrete Comput. Geom. 44 (2010),
        860-882.
    %H  Erdős problem #506, <a href="https://www.erdosproblems.com/506">erdosproblems.com/506</a>
    %e  a(6) = 8: (0,0), (20,0), (5,15), (5,0), (2,6), (10,10) (a triangle and the feet of its altitudes)
        determine 3 four-point circles and 5 three-point circles; the sides are 3-point lines.
    %e  a(8) = 17: (0,0), (0,5), (0,10), (0,15), (3,6), (3,9), (5,5), (5,10) determine 11 four-point circles and
        6 three-point circles.  Segre's projected cube (two concentric squares) determines 18.
    %e  a(9) = 25: (0,0), (25,0), (-25,0), (0,25), (0,-25), (7,24), (-7,-24), (24,7), (-24,-7).
    %Y  Cf. A003035 (orchard problem, maximal number of 3-point lines).
    %K  nonn,hard,more

## 11. Draft comment for teorth/erdosproblems issue #297

> The sequence c(n) of #506 (minimum number of circles through ≥ 3 points of n points not all on a line or a
> circle; collinear triples determine no circle — Elliott's convention, as in the Lean statement) begins
> c(4..10) = 3, 5, 8, 11, 17, 25, 33.  The Elliott–Purdy–Smith value C(n−1,2) + 1 − ⌊(n−1)/2⌋ (= 3, 5, 9, 13, 19,
> 25, 33, …) is an upper bound for all n and is exact for n = 4, 5, 9, 10, but c(6) = 8, c(7) = 11, c(8) = 17
> lie below it; c(8) = 17 also beats Segre's projected cube (18).  Extremal sets: triangle + feet of altitudes
> (n = 6), orthocentric quadrangle + its diagonal points (n = 7), and for n = 8 the integer set (0,0), (0,5),
> (0,10), (0,15), (3,6), (3,9), (5,5), (5,10), obtained by inverting the orthocentric configuration in a point
> common to an altitude, the circumcircle and the nine-point circle.  The lower bounds are computer-assisted
> (enumeration of Möbius block structures under Sylvester–Gallai constraints after inversion, then
> non-realisability of the survivors: Fano, Möbius–Kantor (8_3), an angle lemma, Miquel's theorem, and for
> n = 10 a Gröbner-basis exclusion of two concentric regular pentagons with ten collinear triples); every
> step was reproduced by independent code.  We conjecture c(n) = C(n−1,2) + 1 − ⌊(n−1)/2⌋ for all n ≥ 9
> (searches found nothing better up to n = 16).  Full write-up, certificates and code: `docs/NOTE_506.md`
> and `experiments/506/` in the repository.  This work was carried out by AI agents (Claude) with independent
> re-verification; a human check of the write-up is in progress.

## References

- [Er61] P. Erdős, Some unsolved problems, Magyar Tud. Akad. Mat. Kutató Int. Közl. 6 (1961) 221–254.
- [El67] P. D. T. A. Elliott, On the number of circles determined by n points, Acta Math. Acad. Sci. Hungar.
  18 (1967) 181–188.
- [PS10] G. B. Purdy, J. W. Smith, Lines, circles, planes and spheres, Discrete Comput. Geom. 44 (2010)
  860–882; arXiv:0907.0724.
- [BB94] A. Bálintová, V. Bálint, On the number of circles determined by n points in the Euclidean plane,
  Acta Math. Hungar. 63 (1994) 283–289.
- [Zh11] R. Zhang, On the number of ordinary circles determined by n points, Discrete Comput. Geom. 46 (2011)
  205–211.
- [LMMNS17] A. Lin, M. Makhul, H. N. Mojarrad, E. Naslund, K. Swanepoel, On the number of ordinary circles,
  arXiv:1412.8314.
- [LMMSSZ18] A. Lin, M. Makhul, H. N. Mojarrad, J. Schicho, K. Swanepoel, F. de Zeeuw, On sets defining few
  ordinary circles, Discrete Comput. Geom. 59 (2018) 59–87.
- [KM58] L. M. Kelly, W. O. J. Moser, On the number of ordinary lines determined by n points, Canad. J. Math.
  10 (1958) 210–219.
- [CS93] J. Csima, E. T. Sawyer, There exist 6n/13 ordinary points, Discrete Comput. Geom. 9 (1993) 187–202.
- [GT13] B. Green, T. Tao, On sets defining few ordinary lines, Discrete Comput. Geom. 50 (2013) 409–468.
- [CM68] D. W. Crowe, T. A. McKee, Sylvester's problem on collinear points, Math. Mag. 41 (1968) 30–34.
- [BGS74] S. A. Burr, B. Grünbaum, N. J. A. Sloane, The orchard problem, Geom. Dedicata 2 (1974) 397–424;
  OEIS A003035.
- [Gr09] B. Grünbaum, Configurations of Points and Lines, Grad. Stud. Math. 103, AMS 2009 (the (8_3)
  configuration of Möbius and Kantor is not realisable in the real plane).
- [Mi38] A. Miquel, Théorèmes sur les intersections des cercles et des sphères, J. Math. Pures Appl. 3 (1838)
  517–522.
- google-deepmind/formal-conjectures, `FormalConjectures/ErdosProblems/506.lean`; teorth/erdosproblems,
  issue #297 (May 2026).

## Appendix: file index (repository `Solving-hard-math-problems-`, `experiments/506/`)

- `RESULTS.md` — the table with statuses, certificates and dependencies (this note's companion).
- `circles_exact.py` — exact circle counter (rational arithmetic).
- `theory/REPORT.md`, `theory/*.py` — lemmas, constructions, enumerations for n ≤ 9 (with errata).
- `literature/REPORT.md` — literature survey; `segre_cube_check.py` — exact analysis of the projected cube.
- `verify_independent/` — main agent's independent relaxation (CP-SAT), Angle Lemma and (8_3) certificates,
  complete n = 9 enumeration, n = 10 skeleton and pentagon checks.
- `skeptics/ADJUDICATION.md`, `skeptics/*` — adversarial verification of c(5)…c(9) and the lemmas.
- `n10-continue/REPORT.md`, `n10-continue/*` — the proof of c(10) ≥ 33 (skeletons, two-stage enumeration,
  filters, exact analyses; `runs/post_all.json` holds the 37 classes and the survivor).
- `audit-n10-continue/`, `audit-n10-cpsat/`, `audit-n10-enum/` — independent re-derivations for n = 10.
- `n10-enum/orderly3/` — orderly-generation enumeration (combinatorial bound c(10) ≥ 32).
- `n10-cpsat/v2/` — first CP-SAT attack (partial; its §4.0 hand proof is wrong, see `RESULTS.md` §3.4).
- `search-inversion/`, `search-local/`, `audit-search-*` — construction searches n = 9…16 and their audits.
