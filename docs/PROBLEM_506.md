# Problem dossier — Erdős #506: minimum number of circles determined by n points

## Statement
Let P be a set of n points in the plane, not all on one line and not all on one circle. A circle is
*determined* by P if it passes through at least three points of P (three collinear points determine no
circle). Let c(n) be the minimum number of determined circles over all such P. Determine c(n).

Source: Erdős, "Some unsolved problems" (1961); erdosproblems.com/506; formal statement in
google-deepmind/formal-conjectures `FormalConjectures/ErdosProblems/506.lean` (`erdos_506.variants.small_n`
is marked research-open for 4 ≤ n ≤ 393).

## What was known before this session
- Elliott (1967), corrected by Purdy–Smith (2010): for n > 393,
  c(n) = f(n) := C(n−1,2) + 1 − ⌊(n−1)/2⌋, attained by n−1 points on a circle in antipodal pairs plus the centre.
- Segre: a projected cube shows 8 points can determine fewer than C(7,2) = 21 circles.
- No exact value of c(n) for any 4 ≤ n ≤ 393 was recorded anywhere; no OEIS sequence exists
  (teorth/erdosproblems issue #297, May 2026, asks for exactly this computation).
- Related asymptotic work concerns *ordinary* circles (Bálintová–Bálint 1994, Zhang 2011,
  Lin–Makhul–Mojarrad–Naslund–Swanepoel 2017) and gives nothing quantitative for small n.

## Reformulation used throughout
Work in the Möbius plane (plane plus a point at infinity). A *block* is a circle or line containing ≥ 3
points of P; every triple of P lies in exactly one block; two blocks share ≤ 2 points; two lines share ≤ 1.
Writing D = Σ over blocks of size k ≥ 4 of (C(k,3) − 1) and ℓ = number of lines with ≥ 3 points,

    circles(P) = C(n,3) − D − ℓ.

Inverting about a point p ∈ P turns the blocks through p into the ≥3-point lines of a real (n−1)-point
set, so Sylvester–Gallai (and its quantitative forms) constrains the "derived structure" at every point:
rich blocks through p cover at most C(n−1,2) − o(n−1) pairs, where o(m) ≥ 1 is the minimum number of
ordinary lines of m non-collinear points. Lines cover at most C(n,2) − o(n) pairs.

## Hypotheses and how each is tested
| # | Hypothesis | Test | Status |
|---|-----------|------|--------|
| H1 | The two-concentric-squares configuration (n = 8, 18 circles) beats f(8) = 19 | exact rational count | confirmed; it is Segre's cube example |
| H2 | Exceptions to f(n) exist only at small n (block-size-4 designs cannot beat n²/2 for n ≳ 15) | counting argument; Lemmas A/B/C on the largest block | argued in `experiments/506/theory/REPORT.md` |
| H3 | c(6) = 8, c(7) = 11 (below f) | constructions from orthocentric systems; lower bound = CP-SAT relaxation + Angle Lemma | constructions verified exactly; relaxation reproduced independently; Angle Lemma re-proved (Gröbner) |
| H4 | c(8) = 17 (below f and below Segre's 18) | integer construction; lower bound = CP-SAT relaxation (SG-only caps suffice) | construction verified exactly; relaxation optimum 39 reproduced independently |
| H5 | c(9) = f(9) = 25 | relaxation enumeration + hereditary Sylvester–Gallai at degree-7 points + (8_3) non-realisability | under independent verification |
| H6 | c(n) = f(n) for all n ≥ 9 | n = 10 via CP-SAT with hereditary constraints; structural lemmas for larger n | open |
| H7 | Every 8-point set with ≤ 18 circles carries the cube structure | enumeration of structures with count ≤ 18 | claimed by theory agent, to be re-checked |

## Deliverables
1. Exact values c(4..9) with certificates: explicit integer coordinates for the upper bounds and
   machine-checkable enumerations for the lower bounds (`experiments/506/`).
2. A candidate OEIS sequence (c(n), n ≥ 4): 3, 5, 8, 11, 17, 25, ... and a reply to
   teorth/erdosproblems issue #297.
3. Structural theorems for extremal configurations at n = 8, 9.
