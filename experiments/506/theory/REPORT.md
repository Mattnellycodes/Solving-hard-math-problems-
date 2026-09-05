# Erdős #506 — lower bounds and exact values of c(n) for small n

Directory: `experiments/506/theory/` (all code and logs referenced below live here).
Theory/lower-bound agent, 2026-09-05.  Every claim is labelled **proved**, **computed**
(exact computer verification, code in this directory) or **conjectured**.

## 0. Summary

Let c(n) be the minimum number of circles through ≥ 3 points of an n-point planar set that is not
contained in a line or a circle, and f(n) = C(n−1,2) + 1 − ⌊(n−1)/2⌋ the Elliott/Purdy–Smith value
(proved optimal for n > 393).

| n | c(n) | f(n) | status | extremal / record configuration |
|---|------|------|--------|---------------------------------|
| 4 | 3 | 3 | proved (elementary) | 3 collinear + 1, or 2 chords through the 4th point |
| 5 | 5 | 5 | proved (hand, §3.1; enumeration confirms) | 4 concyclic points, two chords through the 5th |
| 6 | **8** | 9 | proved (§3.2) — **beats the formula** | triangle + feet of its three altitudes |
| 7 | **11** | 13 | proved (§3.3) — **beats the formula** | orthocentric system + its 3 diagonal points |
| 8 | **17** | 19 | proved (§3.4) — **beats the formula, new record** (old: 18) | inversion of the orthocentric configuration in the point where an altitude, the circumcircle and the nine-point circle concur |
| 9 | 25 | 25 | proved (§3.5) | 8 concyclic points + a point on 4 chords (e.g. antipodal pairs + centre) |
| 10 | ? | 33 | open; combinatorial relaxation attacked by CP-SAT (§6) | antipodal configuration gives 33 |

Integer certificates (verified exactly by `certify.py`, an implementation independent of
`../circles_exact.py`, and by `../circles_exact.py`):

* n = 6, 8 circles: `(0,0),(20,0),(5,15),(5,0),(2,6),(10,10)`
* n = 7, 11 circles: the six points above and `(5,5)`
* n = 8, **17 circles**: `(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)`
  (17 = 11 four-point circles + 6 three-point circles; the lines with ≥ 3 points are x = 0 (4 points),
  {(0,0),(3,6),(5,10)} and {(0,15),(3,9),(5,5)}.)

The proofs for n = 6, 7, 8, 9 are computer-assisted case analyses over *abstract Möbius block
structures* (§4), pruned by necessary realisability conditions that are theorems (Sylvester–Gallai
and its hereditary form, non-realisability of the Fano plane and of the (8_3) Möbius–Kantor
configuration — the latter proved from scratch in `eight_three_light.py`), and finished by one
geometric lemma (the Angle Lemma 2.3, with an algebraic certificate).  The enumerations were
cross-validated by three independent programs (orderly generation, a hand-guided enumeration, and a
CP-SAT model).

## 1. Setting: the Möbius reformulation

Identify the plane with the sphere minus a point.  For a finite set P on the sphere (Möbius plane) a
**block** is a circle (of the sphere) containing ≥ 3 points of P; B(P) is the set of blocks, regarded
as subsets of P.  Facts (proved): every triple of P lies in exactly one block; two blocks share ≤ 2
points.  If ∞ ∉ P, the Euclidean lines with ≥ 3 points of P are the blocks through ∞, so

    circles(P) = |B(P)| − ℓ = C(n,3) − D(F) − ℓ,                                       (1)

where F = blocks of size ≥ 4, D(F) = Σ_{B∈F} (C(|B|,3) − 1) (the "deficit"), and ℓ = number of
blocks through ∞ ("lines"); lines pairwise share ≤ 1 point of P.  Hence

    c(n) = min over (F, L) realisable of  C(n,3) − D(F) − |L|,                          (2)

and an equivalent symmetric form: c(n) = min over (n+1)-point Möbius sets Q and q ∈ Q of the number
of blocks of Q avoiding q.  Inverting P in a point p ∈ P maps the blocks through p to straight lines
of the (n−1)-point set P∖{p} ("derived structure at p") and the other blocks to circles.

## 2. Lemmas

**Lemma 2.1 (intersection).** Distinct blocks share ≤ 2 points; distinct lines share ≤ 1 point; a
block of size n means P is concyclic/collinear (excluded).  *Proved (trivial).*

**Lemma 2.2 (derived Sylvester–Gallai, hereditary).** Let o(m) be the minimum number of ordinary
(2-point) lines of m non-collinear points in R²; o(m) ≥ 1 (Sylvester–Gallai), o(3..8) = 3,3,4,3,3,4
(Kelly–Moser for m = 7; Csima–Sawyer o(m) ≥ 6m/13 for m ≠ 7).  For every p ∈ P the blocks of size
≥ 4 through p cover at most C(n−1,2) − o(n−1) pairs of P∖{p}; the lines cover at most
C(n,2) − o(n) pairs of P.  Hereditary form: for every p and every S ⊆ P∖{p} (|S| ≥ 3) not contained
in one block through p, some pair of S lies in no block through p together with a third point of S;
likewise for the lines and every S ⊆ P.  *Proved:* P∖{p} inverted about p is a non-collinear real
point set whose ≥ 3-point lines are exactly the derived blocks; apply Sylvester–Gallai to S.
Consequences used: (i) no derived structure contains a Fano plane on 7 points (a Fano plane covers
all 21 pairs); (ii) at most 6 four-blocks through a point when n = 8 (7 would be a Fano plane).

**Lemma 2.3 (Angle Lemma).** Let L1,…,L4 be four lines in R², pairwise non-parallel, no three
concurrent, with vertices P_ij = L_i ∩ L_j.  Then the three "diagonal quadruples"
{P12,P34,P13,P24}, {P12,P34,P14,P23}, {P13,P24,P14,P23} are not all concyclic.
*Proof.* With directed angles mod π, θ_i the direction of L_i, four points A,B,C,D are concyclic iff
∡(CA,CB) = ∡(DA,DB).  For the first quadruple take C = P12, D = P34, A = P13, B = P24: the condition
is θ2 − θ1 ≡ θ4 − θ3, i.e. θ1 + θ4 ≡ θ2 + θ3.  Likewise the second gives θ1 + θ3 ≡ θ2 + θ4 and the
third θ1 + θ2 ≡ θ3 + θ4 (mod π).  Adding the first two: 2θ1 ≡ 2θ2, so L1 ⊥ L2; adding the first and
third: L1 ⊥ L3; hence L2 ∥ L3, contradiction.  ∎  *Algebraic certificate* (`angle_lemma_check.py`):
in slope coordinates the three concyclicity determinants factor as (degeneracy factors) × A_i with
A1 = m2m3m4 + m2 + m3 − m4, A2 = m2m3m4 + m2 − m3 + m4, A3 = m2m3m4 − m2 + m3 + m4, and
A1 − A2 = 2(m3 − m4); the ideal saturated by the non-degeneracy conditions is the unit ideal.
Consequence: **six points forming a complete quadrilateral (four 3-point lines) cannot have all three
"opposite-pair complements" concyclic** — the unique combinatorial way to reach 7 circles with n = 6.

**Lemma 2.4 ((8_3) is not real).** Every family of 8 triples on 8 points pairwise sharing ≤ 1 point
is isomorphic to the Möbius–Kantor configuration {i,i+1,i+3 mod 8} (*computed*, brute force), and it
has no realisation in RP²: with the frame p0=(1,0,0), p1=(0,1,0), p2=(0,0,1), p5=(1,1,1) (no three of
{0,1,2,5} on a configuration line) the incidences force p3=(1,1,0), p6=(a,1,1), p7=(b,0,1),
p4=(0,c,1) with a = b = 1 − c and c² − c + 1 = 0 (*computed exactly*, `eight_three_light.py`; the
chart choices are exhaustive because the points are distinct).  Hence 8 real points span at most 7
three-point lines, and **no point of a 9-point Möbius set lies in 8 four-blocks**.

**Lemma 2.5 (A: block of size n−1).** If P has a block B of size n − 1 with x ∉ B then
circles(P) ≥ f(n), with equality iff B is a circle and x lies on ⌊(n−1)/2⌋ lines each through two
points of B.  *Proved:* the C(n−1,2) triples {x,a,b} lie in pairwise distinct blocks (a common block
would share 3 points with B); if B is a line no other line exists (it would share 2 points of B and
∞); if B is a circle the lines through x cut B in disjoint pairs, so ℓ ≤ ⌊(n−1)/2⌋.

**Lemma 2.6 (B: block of size n−2).** If the largest block has m = n − 2 points then
circles(P) ≥ m² − m + 1 − 3⌊m/2⌋.  For n ≥ 7 this exceeds f(n) (15 > 13, 22 > 19, 34 > 25, …), so
no configuration beating the formula has a block of size n − 2 when n ≥ 7; for n = 6 the bound is 7
and the truth is 8 (Angle Lemma).  *Proved:* with B, x, y: the only blocks of size ≥ 4 besides B are
{x,y,a,b} with pairwise disjoint pairs {a,b} ⊂ B (t ≤ ⌊m/2⌋ of them); |B(P)| = 1 + 2C(m,2) + m − 3t;
lines: x-lines and y-lines cut B in disjoint pairs and must use "crossing" pairs, plus possibly one
line {x,y,b} which exists only if some b ∈ B is on no x-line; in all cases ℓ ≤ m.

**Lemma 2.7 (C: largest block / chunk lemma).** Let B be a largest block, |B| = m, r = n − m ≥ 1.
Then the number N2 of blocks meeting B in exactly 2 points satisfies
N2 ≥ r·C(m,2) − m·r(r−1)/4, hence
circles(P) ≥ 1 + r·C(m,2) − m·r(r−1)/4 − ℓ,  ℓ ≤ ⌊(C(n,2) − o(n))/3⌋.
*Proof.* For a ≠ b in B the blocks containing {a,b} other than B have exactly two points of B and
partition R = P∖B into "chunks"; let g_ab be their number.  For fixed a, chunks belonging to different
b, b′ share ≤ 1 point (the blocks {a,b}∪S and {a,b′}∪S′ are distinct), so the pairs of R inside chunks
are counted at most once: Σ_{b≠a} Σ_chunks C(|S|,2) ≤ C(r,2).  Since r − g_ab = Σ_chunks (|S| − 1) ≤
Σ_chunks C(|S|,2), summing over a gives 2Σ_{pairs}(r − g_ab) ≤ m·C(r,2), i.e. N2 = Σ g_ab ≥
r·C(m,2) − m·r(r−1)/4.  ∎  (`danger_zone.py` tabulates, for each n ≤ 30, the largest-block sizes not
excluded by Lemmas A, B, C for a configuration with < f(n) circles; e.g. n = 10: m ∈ {4,5,6};
n = 20: m ∈ {4,…,9}; n = 30: m ∈ {4,…,13}.)

**Lemma 2.8 (counting bounds).** (a) Every point lies on ≥ o(n−1) blocks of size exactly 3, so
b₃ ≥ n·o(n−1)/3 and circles(P) ≥ n·o(n−1)/3 + b_{≥4} − ℓ.  (b) Σ_{B∋p} C(|B|−1,2) ≤ C(n−1,2) − o(n−1)
for every p; in particular the number of 4-blocks through p is at most t₃(n−1) (orchard number) and
b₄ ≤ n·t₃(n−1)/4.  (c) ℓ ≤ ⌊(C(n,2) − o(n))/3⌋, and ℓ = number of ≥ 3-point lines ≤ t_{≥3}(n).
*Proved* (immediate from Lemma 2.2).

**Lemma 2.9 (Miquel closure; stated, not needed).** If 8 points carry 5 of the 6 "faces" of a
combinatorial cube as concyclic quadruples, the 6th is concyclic (Miquel's theorem in the real
inversive plane).  Any *exact* structure must be Miquel-closed; this can be used as an extra pruning
rule (it was not needed for n ≤ 9).

## 3. Exact values

### 3.1 n = 4, 5 (proved by hand)
n = 4: circles = 4 − (collinear triples), at most one collinear triple, so c(4) = 3.
n = 5 (C(5,3) = 10): a 4-block B and x: 6 further blocks of size 3; if B is a line, ℓ = 1 and 6
circles; if B is a circle, ℓ ≤ 2 (disjoint chords through x) and circles ≥ 5, attained.  No 4-block:
10 blocks, ℓ ≤ 2 (partial STS(5)), ≥ 8.  So c(5) = 5.  (Enumeration with target 4: no structure.)

### 3.2 n = 6: c(6) = 8 (proved)
Construction: triangle A,B,C and the feet D,E,F of its altitudes: three 3-point lines (sides), three
4-point circles (circles with diameters AB, BC, CA each contain two feet), and 5 ordinary circles:
20 = 3 + 12 + 5 triples, 8 circles (`certify.py`).
Lower bound: by (1), circles = 20 − 3q − ℓ where q ≤ 3 is the number of 4-blocks (no 5-block: Lemma A
gives ≥ 9).  q ≤ 2 gives ≥ 20 − 6 − 4 = 10 (ℓ ≤ 4 for 6 points).  q = 3 forces the three 4-blocks
to be the complements of a perfect matching {12|34|56}; the 8 uncovered triples are the transversals;
ℓ = 4 would need four transversal lines (a 4-block line kills all transversals), i.e. a complete
quadrilateral whose three opposite-pair complements are concyclic — impossible by the Angle Lemma.
So circles ≥ 20 − 9 − 3 = 8.  (The enumerator with target 7 finds exactly this one candidate.)

### 3.3 n = 7: c(7) = 11 (proved)
Construction: orthocentric system A,B,C,H (H the orthocentre) with the three diagonal points D,E,F of
the complete quadrangle: 6 lines (the sides of the quadrangle), 6 four-point circles (for each pair
X,Y ∈ {A,B,C,H} the circle with diameter XY contains two of D,E,F), 5 ordinary circles;
35 = 6 + 24 + 5; 11 circles (`certify.py`).
Lower bound: circles ≤ 10 needs D + ℓ ≥ 25 with ℓ ≤ 6 (7 lines on 7 points would be a Fano plane,
excluded by Sylvester–Gallai), so D ≥ 19; blocks of size ≥ 5 give ≥ 13 (Lemmas A, B); hence 7
four-blocks pairwise sharing ≤ 2 points, whose complements are 7 triples pairwise sharing ≤ 1 point:
the Fano plane.  In this family every point p lies in 4 blocks and avoids 3; inverting about p gives
six points with 4 three-point lines (a complete quadrilateral) and 3 four-point circles which are
exactly the three opposite-pair complements — impossible by the Angle Lemma.  So c(7) ≥ 11.
(Enumeration with target 10: exactly this one candidate, combinatorial count 8.)

### 3.4 n = 8: c(8) = 17 (proved; new record configuration)
**Construction (proved/computed).** Take A = (0,3), B = (1,0), C = (3,0); orthocentre H = (0,−1),
feet D = (6/5,−3/5), E = (2,1), F = (0,0).  The point O = (0,1) is simultaneously the reflection of H
in BC (hence on the circumcircle of ABC) and the midpoint of AH (hence on the nine-point circle), and
lies on the altitude line AHF.  The 8-point Möbius set {A,B,C,H,D,E,F,∞} has the cube structure
(12 four-blocks, 8 three-blocks; the 6 lines are the 6 four-blocks through ∞); O lies on three of
its blocks (in cube language: a diagonal plane through ∞ and the two corner triples N(v), N(w)
of an edge vw lying on that plane — the 48 combinatorial line-triples of the cube structure come in
three kinds, computed: 24 "face + corner triples of an edge of the opposite face", 12 "diagonal plane
+ two adjacent corner triples on it", 12 "diagonal plane + two antipodal corner triples on it"), so
inverting in O yields 8 finite points with 20 − 3 = 17 circles:
(0,3/2),(1/2,1/2),(3/10,9/10),(0,1/2),(3/10,3/5),(1/2,1),(0,0),(0,1), i.e. after scaling by 10
the integer set (0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10).  Verified exactly by two
independent programs and the float code of `../explore_families.py`.  The whole family: triangles
with A = (0,3), B = (b,0), C = (3/b,0), 0 < b < √3 (up to similarity, a 1-parameter family), all give
17 generically (*computed* for b = 1; conjectured for all b).
**Lower bound (computer-assisted, `mobius_enum.py 8 17`, 166 nodes, < 1 min).** The enumeration of
all abstract structures (F, L) satisfying Lemmas 2.1–2.2 with C(8,3) − D − ℓ ≤ 17 returns exactly one
structure: the 12-block cube structure (SQS(8) minus a parallel class; all degrees 6) with ℓ_max = 3
and count exactly 17.  Since no structure has count ≤ 16, c(8) ≥ 17.  The same conclusion holds
using Sylvester–Gallai alone (o(m) ≥ 1) instead of Kelly–Moser/Csima–Sawyer (validation run V4).
Independent confirmations: (V2) `verify_n8_independent.py` — starting from a point of degree 6, whose
derived structure must be "Fano minus a line" (unique), all 11- and 12-block families are built
directly: 12-block: 1 family (the cube), ℓ_max = 3; 11-block: 13 labelled families, ℓ_max = 4 (count
19); families with a block of size ≥ 5 never reach 17; (V5) the CP-SAT model `cpsat_model.py 8`
proves the combinatorial optimum D + ℓ = 39, i.e. 17, in 24 s.
By-product (structure theorem, *computed*): the only abstract structure with count ≤ 18 besides the
cube is a 10-block structure needing 8 three-point lines on 8 points (Möbius–Kantor, impossible by
Lemma 2.4).  Hence **every 8-point set with ≤ 18 circles carries the cube structure** (12 four-point
blocks), and the sets with 17 circles are exactly the cube-structure sets with ∞ on three blocks.

### 3.5 n = 9: c(9) = 25 (proved, computer-assisted)
Upper bound: 8 points on a circle in 4 antipodal pairs + the centre (25, `certify.py`).
Lower bound: circles ≤ 24 needs D + ℓ ≥ 60, ℓ ≤ 10, D ≥ 50.  Blocks of size 8, 7, 6 are excluded by
Lemmas A, B, C (≥ 25, 34, 28).  *Case b₅ = 0:* D = 3b₄ ≥ 50 forces b₄ ≥ 17 and a point in ≥ 8
four-blocks, whose derived structure is an (8_3) — impossible by Lemma 2.4.  *Case b₅ ≥ 1:* the
two-phase enumerator `mobius_enum2.py 9 24 --big 5,6,7,8` (4981 nodes, 5 s) returns exactly one
structure: two 5-blocks sharing a point plus 12 four-blocks (degrees 2,7,…,7), ℓ_max = 6, count 24.
`hereditary_sg.py` shows that at every point of degree 7 the derived structure contains a Fano plane
on the other 7 points (the four points of one 5-block form a complete quadrangle whose diagonal
points would be collinear), contradicting Lemma 2.2.  Hence c(9) ≥ 25.
Structure theorem (*computed*, target 25): the only abstract structures with count ≤ 25 are the
excluded one above and the antipodal structure (one 8-block, 4 lines).  So **every 9-point set with
25 circles consists of 8 concyclic points and a point lying on 4 chords** (the chords need not be
diameters; the 9th point need not be the centre).

## 4. Enumeration methodology and validation

`mobius_enum.py` (n ≤ 8): Read–Faradžev orderly generation of families F of blocks of size 4..n−1
pairwise sharing ≤ 2 points; canonical = lexicographically minimal sorted block list under S_n
(vectorised permutation table).  Appendix A proves that the parent (remove the largest block) of a
canonical family is canonical, so each isomorphism class is visited once.  Pruning: (C2) coverage
caps, and D(F) + potential ≥ D_min where potential ≤ (remaining coverage budget)/3 and ≤ Σ deficits of
compatible candidates.  For each F, `line_sets` enumerates maximal sets of pairwise ≤ 1-intersecting
blocks (from F and the uncovered triples) within the pair budget (C3).
`mobius_enum2.py` (n = 9, 10): phase 1 enumerates canonical big-block families (sizes ≥ 5), phase 2
extends by 4-blocks with canonicity tested only against Aut(F_big) (valid because the canonical form
compares big blocks first; Appendix A).
Validation: (V1) `validate_canon.py` — orderly generation reproduces brute-force isomorphism-class
counts for n = 5, 6, 7 (2, 5, 18 classes); (V2) independent hand-guided enumeration for n = 8;
(V3) independent float recount of the 17-configuration; (V4) n = 8 with Sylvester–Gallai-only caps
(same unique candidate); (V5) `cpsat_model.py`: CP-SAT optimum of D + ℓ under (C1)–(C3) is 13, 27,
39 for n = 6, 7, 8 (combinatorial minima 7, 8, 17), matching the enumerations.
`cpsat_hereditary.py` adds the hereditary constraints (no derived Fano plane, no derived (8_3),
also for the lines) as linear constraints and is the tool proposed for n = 10 (§6).
Logs/JSON: `enum_n8_t17.json`, `enum_n8_t18.json`, `enum2_n9_t24.json`, `enum2_n9_t25.json`,
`cpsat_n9.log`, `cpsat_n10.log`.

## 5. What a configuration beating the formula must look like (general n)

All *proved* unless marked.
1. Its largest block has size m with 4 ≤ m ≤ m*(n), where m*(n) is given by Lemmas A, B, C
   (`danger_zone.py`): m ≤ n − 3 for n ≥ 7, and roughly m ≲ n/2 − O(1) (table: n = 10 → m ≤ 6,
   n = 15 → m ≤ 8, n = 20 → m ≤ 9, n = 30 → m ≤ 13).  Nothing with a big block can beat f(n).
2. It needs D + ℓ ≥ C(n−1,3) + ⌊(n−1)/2⌋, with ℓ ≤ ⌊(C(n,2) − o(n))/3⌋ ≈ n²/6 and, by Lemma 2.8(b),
   D ≤ Σ_B C(|B|,3) − b_{≥4} ≤ n(C(n−1,2) − o(n−1))/3 − b_{≥4}.  In words: **almost every point must
   lie on the maximum possible number of rich derived lines**, i.e. every derived structure must be a
   near-extremal orchard configuration, and the line structure at ∞ must be near-extremal too.  For
   n ≤ 9 this is exactly what kills the candidates (Fano / (8_3) sub-configurations).
3. It has no block of size ≥ 5 for n ≤ 9 (computed), and every point lies on ≥ o(n−1) ordinary
   circles (3-blocks).
4. Heuristic (conjectured): for n ≥ 10 the deficit obtainable from 4- and 5-blocks is at most
   ≈ 3n·t₃(n−1)/4 + O(n), which falls below C(n−1,3) − n²/6 for n ≥ 10, so we **conjecture c(n) =
   f(n) for all n ≥ 9**, and that n = 6, 7, 8 are the only exceptions.  The CP-SAT relaxation for
   n = 10 (§6) is the natural first test.

## 6. n = 10 and the plan

Reduction (proved): circles ≤ 32 needs D + ℓ ≥ 88, ℓ ≤ 13 (Csima–Sawyer o(10) ≥ 5), D ≥ 75;
largest block ≥ 7 excluded (Lemma C gives ≥ 41, Lemma B 45, Lemma A 33); two 6-blocks impossible
(they would cover all 10 points and admit ≤ 12 four-blocks).  Remaining: largest block 4, 5 or 6.
* Largest block 4: b₄ ≥ 25 forces every point in exactly 10 four-blocks and ℓ ≥ 13, i.e. 13
  three-point lines on 10 points; this is impossible by the orchard value t₃(10) = 12 (Burr–Grünbaum–
  Sloane 1974, *cited*), or, self-contained, by enumerating the 13-triple systems on 10 points and
  checking realisability as in Lemma 2.4 (planned).
* Largest block 5 or 6: `cpsat_hereditary.py 10` (CP-SAT with (C1)–(C3) + hereditary Fano/(8_3)
  constraints) — see `cpsat_n10.log` and §7 for its status.  If the optimum is D + ℓ ≤ 87 the value
  c(10) = 33 follows from theorems only; if a structure with 88 or more appears, it is to be attacked
  as in §3.5 (hereditary SG on subsets, Miquel closure, then exact realisation equations via Gröbner
  bases with a Möbius frame of 3 points).

Plan to certify any small n (implemented for n ≤ 9, scalable to n ≈ 11–12):
1. Combinatorial relaxation: enumerate all (F, L) with C(n,3) − D − ℓ ≤ target under Lemmas 2.1–2.2
   (orderly generation or CP-SAT); use the largest-block lemmas to restrict block sizes.
2. Hereditary pruning: derived structures at every point and the line structure at ∞ must satisfy
   hereditary Sylvester–Gallai (`hereditary_sg.py`), have ≤ t₃(m) three-point lines, and be Miquel-
   closed.
3. For each survivor: try to realise numerically (random starts + Newton on the incidence equations,
   3 points fixed by a Möbius transformation); a numerical solution is then certified exactly (rational
   or algebraic coordinates, `certify.py`).  Otherwise prove non-realisability exactly: incidence
   polynomials with the non-degeneracy conditions saturated (Rabinowitsch), Gröbner basis = (1) over
   C, or a univariate eliminant with no real roots (as for (8_3)), or a geometric argument (Angle
   Lemma type).  Realisations of "at least S" reduce to structures with more incidences, which are
   in the list, so this is sound for lower bounds.
Cost: n = 8 needs seconds, n = 9 seconds (two-phase), n = 10 minutes–hours (CP-SAT), n = 11–12
would need symmetry-broken CP-SAT or a local-structure enumeration (derived structures must be
extremal orchard configurations, which are classified for ≤ 12 points).

## 7. Negative results and dead ends

* The regular cube (two concentric squares) cannot give 17: for none of the three kinds of
  combinatorial line-triples are the three circles concurrent on the sphere (antipodal corner triples
  lie in parallel planes; for the other two kinds the intersection line of the two corner planes
  misses the third circle) — checked by hand; only special members of the 2-parameter cube family
  (orthocentric systems with A = (0,3), B = (b,0), C = (3/b,0)) have three concurrent blocks.
* Adding a 9th point to a cube-structure set is bad (≥ 33 circles): the cube's blocks are already
  Miquel-closed and a new point lies on ≤ 2 of them.
* For n = 9 nothing beats the antipodal configuration; the only combinatorial competitor (two 5-blocks
  through a point + 12 four-blocks) is a "Fano-derived" structure.
* For n = 8 the 11-four-block structures reach only ℓ = 4 (count 19), never 17 or 18.
* Structures with blocks of size ≥ 5 are never competitive for n ≤ 9.
* The web-search budget was exhausted before literature checks; Kelly–Moser, Csima–Sawyer,
  Burr–Grünbaum–Sloane values are quoted from memory; only o(7) = 3, o(8) ≥ 4 and the (8_3)/Fano
  facts are load-bearing, the (8_3) and Fano facts are proved here from scratch, and the n = 8
  result was re-derived with Sylvester–Gallai alone (V4).

## 8. Files

* `mobius_enum.py`, `mobius_enum2.py` — enumerators; `validate_canon.py`, `verify_n8_independent.py`
  — validation; `cpsat_model.py`, `cpsat_hereditary.py` — CP-SAT formulations; `hereditary_sg.py`
  — hereditary SG checker; `angle_lemma_check.py`, `eight_three_light.py` — algebraic certificates;
  `danger_zone.py` — Lemma A/B/C table; `certify.py` — independent exact certificate checker with
  the record coordinates; `enum*.json`, `*.log` — outputs.

## Appendix A — orderly generation

Canonical form: the sorted list of block codes (bitmasks), lexicographically minimal over S_n
(for the two-phase version: the pair (sorted big blocks, sorted 4-blocks) compared big blocks first).
Claim: if F is canonical then F′ = F ∖ {max block} is canonical.  Proof: if σ(F′) < F′ (sorted
lists), the first index i where they differ has σ(F′)_i < F′_i; inserting σ(b) into σ(F′) and
b = max F into F′ (which goes to the end) gives σ(F)_j ≤ σ(F′)_j = F_j for j < i and σ(F)_i ≤
σ(F′)_i < F_i, so σ(F) < F, contradicting canonicity.  In the two-phase version, permutations that
move F_big to a lex-larger family are irrelevant and those making it smaller do not exist, so only
Aut(F_big) matters.  The DFS adds blocks in increasing code order, prunes non-canonical children, and
is therefore exhaustive over isomorphism classes; validated against brute force (V1).

## Addendum — status of the CP-SAT runs and extra checks (end of session)

* `cpsat_model.py 6/7/8` (plain relaxation, (C1)–(C3) only): OPTIMAL, D + ℓ = 13, 27, 39, i.e.
  combinatorial minima 7, 8, 17 — matching the orderly enumeration (V5).
* `cpsat_model.py 10 600` (plain relaxation, 10 min, 2 workers, machine load ≈ 15): FEASIBLE only;
  best structure found = the antipodal one (one 9-block + 4 lines, D + ℓ = 87, 33 circles); dual
  bound still 112, so **inconclusive** as a lower bound for n = 10 (`cpsat_n10_plain.log`).
* `cpsat_hereditary.py 9/10` (with derived-Fano/(8_3) constraints): models built (18360 / 14400
  extra constraints) but neither solve finished within the CPU budget under load; stopped.  Not
  needed for n = 9 (settled by the two-phase enumeration + Lemmas 2.2, 2.4), so the n = 9 result does
  not depend on them.
* `robustness_n9.py`: the n = 9 enumeration re-run with Csima–Sawyer-only bounds (cap_derived 24,
  cap_lines 31) and with Sylvester–Gallai-only bounds (27, 35) returns the same single candidate
  (two 5-blocks through a point + 12 four-blocks), killed by hereditary SG at every degree-7 point.
* `family17.py`: b = 1, 1/2, 3/2, 2/3, 5/4, 1/3, 7/5 all give 8-point sets with exactly 17 circles
  (one 4-point line, two 3-point lines, eleven 4-point circles, six 3-point circles).
* Line-triples of the cube structure: 48, in three kinds (24 face-type, 12 diagonal-adjacent,
  12 diagonal-antipodal); the record configuration realises the diagonal-adjacent kind (the altitude
  line AHF∞, the circumcircle ABC and the nine-point circle DEF; ∞ = (1,1,1) in cube coordinates,
  bipartition {A,B,C,H} | {D,E,F,∞}).
