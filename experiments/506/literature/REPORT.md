# Erdős #506 — literature scout report (2026-09-05)

Scope: everything findable about the exact minimum number c(n) of circles determined by n points
(not all on one line, not all on one circle) for *small* n (n ≤ 393), plus the neighbouring
"ordinary circle" / "circle orchard" literature. Access limits: only search-engine snippets and
github.com pages were reachable; the papers of Elliott (1967) and Bálintová–Bálint (1994) are
paywalled and were NOT read in full — statements about their contents are marked (S) = from
snippets / citing papers, (V) = verified by us, (U) = unverified.

Files in this directory:
- `506.lean` — verbatim copy of `FormalConjectures/ErdosProblems/506.lean` (google-deepmind/formal-conjectures, main).
- `segre_cube_check.py`, `segre_cube_check.out` — exact (sympy) analysis of the cube on its circumsphere.
- `segre_cube_float.out` — float checks of several stereographic projections (uses `evaluate` from `../explore_families.py`).
- `two_squares_exact.out` — exact rational counts for the two-concentric-squares family.

---

## 1. Status of the problem (what the sources actually say)

| Source | Statement |
|---|---|
| erdosproblems.com/506 (S) | "What is the minimum number of circles determined by any n points in R², not all on a circle?" Notes a non-degeneracy condition is intended (not all on a line, or no three on a line). "Resolved by Elliott, who proved that (assuming not all points are on a circle or a line), provided n > 393, the points determine at least C(n−1,2) distinct circles. The problem appears to remain open for small n. Segre observed that projecting a cube onto a plane shows that the lower bound C(n−1,2) is false for n = 8." Problem dated 1961 (Erdős, *Some unsolved problems*, Magyar Tud. Akad. Mat. Kutató Int. Közl. 6 (1961) 221–254). |
| erdosproblems.com status (S, via formal-conjectures issue #4952, auto-sync) | The site now lists 506 as **solved** (repo annotation still "open" → mismatch issue opened by the bot). "Solved" evidently refers to Elliott's large-n theorem; the site text itself still says small n "appears to remain open". Note the site quotes the *uncorrected* bound C(n−1,2). |
| erdosproblems.com/forum/thread/506 (S) | Exists; the snippet shows only the problem text — no additional mathematical content surfaced (could not be fetched). |
| Lean file 506.lean (V, copied here) | `erdos_506` (all n ≥ 4): research open. `variants.large_n` (n > 393): solved, answer C(n−1,2) + 1 − ⌊(n−1)/2⌋, "Elliott [El67], with the correction of Purdy and Smith [PuSm], also reported in [BaBa94]". `variants.small_n` (4 ≤ n ≤ 393): research open. `variants.segre_n_eq_eight`: ∃ 8 points with < C(7,2) = 21 circles. Non-degeneracy used: not collinear and not cospherical (Elliott's weaker condition). Reference list says "[PuSm] Purdy and Smith. No reference found." — see §3, the reference does exist. |
| formal-conjectures (V) | Issue #2128 (AI-generated, closed) → PR #4396 by theebayuser, merged 2026-08-14 (reviewer mo271 asked to separate the solved large-n case from the open small-n case). Issue #4952 (open, bot): status mismatch open/solved. |
| teorth/erdosproblems (V) | **Issue #297 (2026-05-14, label "help wanted", no comments as of today): "Computing sequence for Erdős problem #506 (min circles determined by n points)"** — states that an OEIS search found *no* sequence for f(3), f(4), … and asks for a human-written computation before OEIS submission (parallel to #89/A186704 and #669/A003035). No PR in that repo mentions 506. The "AI contributions" wiki has no entry for 506. |

**Bottom line:** no source gives any exact value of c(n) for any n ≤ 393 beyond "c(8) < 21". No tabulation exists anywhere we could find, and there is no OEIS sequence (confirmed by issue #297 and by our own searches). c(n) for 4 ≤ n ≤ 393 is genuinely open in the literature.

## 2. Elliott 1967

P. D. T. A. Elliott, *On the number of circles determined by n points*, Acta Math. Acad. Sci. Hungar. 18 (1967) 181–188 (paywalled).

What it proves (S, consistently reported by Purdy–Smith 2010, Zhang 2011, Lin et al. 2018, Lin–Swanepoel 2019):
- n points, not all on a line or circle, determine at least C(n−1,2) circles **for n > 393** (i.e. n ≥ 394). Purdy–Smith: "the result is slightly wrong; the correct lower bound is 1 + C(n−1,2) − ⌊(n−1)/2⌋ … Elliott's proof can be modified to show the correct result with a lower bound of 394 for n."
- At least 2n²/63 − O(n) **ordinary** circles (circles through exactly three points). Elliott "introduced the question of ordinary circles" and improved Motzkin (1951, *The lines and planes connecting the points of a finite set*, Trans. AMS 70) who had earlier considered circles determined by point sets.
- Elliott himself cited Segre for an eight-point counterexample (the cube) to the bound C(n−1,2) (S, Purdy–Smith).

Method (S/U): the proof works by inversion — inverting about a point p ∈ P turns the circles through p into lines of the inverted set of n−1 points, so the Kelly–Moser theorem (≥ 3n/7 ordinary lines, 1958) yields ordinary circles through p, and a counting argument over all p gives both the ordinary-circle bound and the total bound. **Why exactly 393 is not documented anywhere accessible**; it is presumably where the explicit linear error terms of this inversion/Kelly–Moser counting are beaten by the main term. We did not see the paper, so the structure of the argument (in particular whether the 393 is an artefact of the 3/7 constant and could be lowered using Csima–Sawyer 6/13 or Green–Tao n/2) is an open question worth checking by someone with journal access. Bálintová–Bálint 1994 and Zhang 2011 improve the ordinary-circle constant by refining exactly this kind of argument, so a re-run of Elliott's total-count argument with modern inputs would very likely lower the 393 threshold — but nobody seems to have done it.

## 3. The Purdy–Smith correction — it IS published

George B. Purdy and Justin W. Smith, *Lines, Circles, Planes and Spheres*, Discrete Comput. Geom. 44 (2010) 860–882; arXiv:0907.0724 (2009). The Lean file's "[PuSm] No reference found" is just a missing citation; this is the reference (the phrase "unpublished corrections" is wrong).

Content relevant here (S):
- Corrected bound: n points not all on a line/circle determine ≥ 1 + C(n−1,2) − ⌊(n−1)/2⌋ ≥ 1 + (n−1)(n−3)/2 circles for n ≥ 394. Extremal example: n−1 points on a circle plus a point p off it; when p is the centre and the n−1 points come in antipodal pairs, the ⌊(n−1)/2⌋ triples (p, x, −x) are collinear and determine no circle. (This is exactly why the "not all on a line" rather than "no three collinear" condition matters: under "no three collinear" the answer would revert to C(n−1,2) for large n.)
- Sphere analogue: n points in R³ not all on a plane/sphere determine ≥ 1 + C(n−1,3) − t₃^orchard(n−1) spheres, best possible under the hypothesis (t₃^orchard = max number of 3-point lines).
- They also give new lower bounds for lines and circles (Beck-type, for point sets with no rich line).

## 4. Bálintová–Bálint 1994

A. Bálintová, V. Bálint, *On the number of circles determined by n points in the Euclidean plane*, Acta Math. Hungar. 63 (1994) 283–289 (paywalled).
- (S) Improves Elliott's ordinary-circle bound from 2n²/63 − O(n) to 11n²/247 − O(n) (= (22/247)·C(n,2) − O(n)).
- (U) The Lean file says the corrected count C(n−1,2)+1−⌊(n−1)/2⌋ is "also reported in [BaBa94]". We could not confirm this from any snippet. If true, the correction predates Purdy–Smith by 16 years. Worth checking by someone with access.

## 5. Segre's cube example — exact analysis (V)

Sources only say "projecting a cube onto a plane" (erdosproblems.com; Elliott citing Segre). The right reading: the 8 vertices (±1,±1,±1) lie on the sphere x²+y²+z²=3; circles on the sphere are plane sections; **stereographic** projection from a pole N on the sphere maps circles through N to lines and all other circles to circles. So the circle count of the projected set is (#planes through ≥3 vertices) − (#such planes through N).

Computed exactly with sympy (`segre_cube_check.out`):
- Planes through ≥ 3 vertices: **20** = 12 four-point planes (6 faces + 6 "diagonal" planes x=±y, y=±z, x=±z, each through two opposite edges) + 8 three-point planes (x±y±z = ±1, each through the three neighbours of a vertex).
- Every non-vertex point of the sphere lies on **at most 2** of these planes (54 such degree-2 points, all listed in the .out file; e.g. (0,0,±√3) on the 4-fold axis, (−1/3,−1/3,5/3), ((−1+√5)/2, (−1−√5)/2, 0)). There is no point of degree ≥ 3.
- Hence the projections give exactly: generic pole → 20 circles; pole on one face circle → 19; pole on any degree-2 point → **18**, and 18 is the minimum over the whole Möbius class of the cube.

Closed-form coordinates of the best projection (pole (0,0,√3), image plane z=0):
```
outer square: (±(3+√3)/2, ±(3+√3)/2)     inner square: (±(3−√3)/2, ±(3−√3)/2)
```
i.e. two concentric, equally oriented squares with the two diagonals y=±x each carrying 4 points; ratio of half-sides = 2+√3. Float check (`segre_cube_float.out`): euclid_count = 18, nblocks = 20.

Relation to our record: our {(±1,0),(0,±1),(±2,0),(0,±2)} is this configuration rotated by 45° with ratio 2 instead of 2+√3. Exact rational checks (`two_squares_exact.out`) give 18 circles (10 four-point circles + 8 ordinary circles, 8 collinear triples on 2 lines) for every ratio tried (2, 3, 5, 7, 5/2, 7/3, 10/3, 37/10). Geometrically, ratio r corresponds to a square prism a×a×c inscribed in a sphere (same 12 four-point planes as the cube); the cube is the special case r = 2+√3. So **our 18-circle configuration is Segre's example** (same combinatorial type, Möbius-equivalent to a square prism on a sphere, projected from the 4-fold axis) — the projection along the 4-fold axis is precisely a pole of maximal degree.

Consequences:
- c(8) ≤ 18 < 19 = C(7,2)+1−⌊7/2⌋: Segre's cube beats not only the uncorrected bound 21 but also the *corrected* formula, so the Purdy–Smith formula is not the true minimum at n = 8. No source states the number 18 explicitly; erdosproblems.com and the Lean file only assert "< 21". Whether 18 is optimal for n = 8 is not addressed anywhere.
- The cube cannot give fewer than 18 (max pole degree 2, exact).

## 6. Related results that bound c(n) from below or constrain extremal sets

1. **Ordinary circles** (a circle through exactly 3 points; all such circles are distinct, so c(n) ≥ #ordinary circles):
   - Elliott 1967: ≥ 2n²/63 − O(n). Bálintová–Bálint 1994: ≥ 11n²/247 − O(n). Zhang 2011 (*On the number of ordinary circles determined by n points*, DCG 46 (2011) 205–211): improves 22/247·C(n,2) to (1/9)·C(n,2) − O(n).
   - Lin, Makhul, Mojarrad, Naslund, Swanepoel, *On the number of ordinary circles* (arXiv:1412.8314, 2014/2017): ≥ n²/4 − O(n) via Green–Tao; n²/4 best possible for even n.
   - Lin, Makhul, Mojarrad, Schicho, Swanepoel, de Zeeuw, *On sets defining few ordinary circles*, DCG 59 (2018) 59–87, arXiv:1607.06597: exact minimum number of ordinary circles for all sufficiently large n: n²/4 − 3n/2 (n ≡ 0 mod 4), n²/4 − 3n/4 + 1/2 (1 mod 4), n²/4 − n (2 mod 4), n²/4 − 5n/4 + 3/2 (3 mod 4). **Structure theorem:** if P spans ≤ Kn² ordinary circles (n large in terms of K) then all but O(K) points lie on an algebraic curve of degree ≤ 4 (lines, circles, ellipses, circular cubics, and the relevant quartics — images of lines/conics/cubics under inversion). Also the **circle orchard problem**: max number of 4-point circles = n³/24 − n²/4 + O(n), exact for large n depending on n mod 8 (e.g. n³/24 − n²/4 + 5n/6 for n ≡ 0 mod 8), attained by circles only.
   - Lin–Swanepoel, *Ordinary hyperspheres and spherical curves* (arXiv:1905.09639, Adv. Geom. 2021): higher-dimensional analogues.
   - All of these are asymptotic ("n sufficiently large") and say nothing quantitative for n ≤ 393. The structure theorem does say, heuristically, that any configuration with O(n²) blocks is "mostly on a degree-≤4 curve", supporting the heuristic that exceptions to the formula live at small n.

2. **Inversion + Sylvester–Gallai-type inputs** (the tool behind Elliott's proof, usable at every n):
   For p ∈ P with P∖{p} not all on one circle/line, invert about p: blocks (circles or lines) through p ↔ lines of an (n−1)-point set not all collinear. Hence, for every such p,
   - deg(p) ≥ n−1 (de Bruijn–Erdős), and
   - the number of size-3 blocks through p ≥ (number of ordinary lines of an (n−1)-set) ≥ 3(n−1)/7 (Kelly–Moser 1958), ≥ 6(n−1)/13 for n−1 ≠ 7 (Csima–Sawyer 1993), ≥ (n−1)/2 for n−1 large (Green–Tao 2013); exact small values of the ordinary-line minimum are known for all small m (Crowe–McKee-type tables).
   If P∖{p} *is* concyclic (p is the one point off the circle C), then all blocks through p are 3-point blocks except that pairs of points of C collinear with p merge — this is the extremal situation of the formula.
   Elliott's total-count bound is a bookkeeping of these facts; re-deriving it with the sharper modern ordinary-line bounds is the obvious route to lowering 393 and is, as far as we can tell, not in the literature.

3. **Combinatorial (design-theoretic) bounds are much weaker.** Blocks form a "3-linear space" (every triple in exactly one block). The hypergraph de Bruijn–Erdős theorem (Alon–Mellinger–Mubayi–Verstraëte, *The de Bruijn–Erdős theorem for hypergraphs*, Des. Codes Cryptogr. 2012, arXiv:1007.4150) gives only ≳ n^{3/2} blocks for such 3-covers, with finite inversive (Möbius) planes of order q (q²+1 points, q³+q circles of size q+1) as the near-extremal examples. Real point sets cannot approach this because every derived structure (blocks through p) must be a real line arrangement — e.g. SQS(8) is unrealisable since its derived design is the Fano plane. So the quadratic lower bound is a genuinely geometric fact; any small-n lower-bound proof must use Sylvester–Gallai-type geometry on derived structures.

4. **Neighbouring Erdős problems** (same tag): #104 (unit circles through ≥3 points, o(n²)?), #827 / #831 (triangles with distinct circumradii), #669 (orchard, OEIS A003035), #89 (distinct distances, A186704). Also Erdős–Purdy questions on the maximum number of 4-point circles (now settled asymptotically by the circle-orchard theorem above).

## 7. Formula values vs. what is known for small n

c_formula(n) = C(n−1,2)+1−⌊(n−1)/2⌋ (proved minimal only for n ≥ 394):

| n | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|---|
| formula | 3 | 5 | 9 | 13 | 19 | 25 | 33 | 41 | 51 |
| literature | — | — | — | — | < 21 (Segre) | — | — | — | — |
| ours (V) | | | | | 18 = Segre | | | | |

Trivial exact value (V, elementary): c(4) = 3 — four points not all concyclic/collinear have at most one collinear triple (two collinear triples share two points and force all four collinear), and the remaining ≥ 3 triples give distinct circles (two triples on one circle would make all four concyclic). So c(4) = 3 = formula. Nothing beyond this is proved in the literature for n ≤ 393.

## 8. Negative results / dead ends (recorded)
- No OEIS sequence for c(n) (own searches + teorth/erdosproblems issue #297). No table of c(n) for any small n anywhere.
- No accessible text explains the origin of 393. Elliott 1967 and Bálintová–Bálint 1994 could not be read (paywall; non-GitHub fetches blocked).
- Erdős's 1961 problem paper (renyi.hu/~p_erdos/1961-22.pdf) not fetchable; whether Segre's remark appears there or only in Elliott's paper is unresolved.
- The Purdy–Smith arXiv PDF (0907.0724) is not fetchable here (proxy 403); everything about it is from snippets.
- Erdősproblems forum thread 506 has no visible mathematical content in snippets.
- Web search budget for this session was exhausted (200 queries) partway; the queries on "n ≥ 394", Bálintová–Bálint abstract, Erdős 1961 text, /history/506 and Sylvester–Gallai-for-circles did not run.

## 9. Bug found in `experiments/506/explore_families.py` (important)
`evaluate()` filters candidate inversion centres with an *exact* 1e-7 bucket match against the points of P (`(key(x),key(y)) in ptkeys`), but the intersection points are computed from block parameters already rounded to 1e-7, so a point of P frequently fails the match and is accepted as a centre. Example (segre_cube_float.out): for the cube image it returns best_count = 11 with best_O = (−2.366, −2.366), which is the point (−(3+√3)/2, −(3+√3)/2) ∈ P (11 = 20 − 9 blocks through a vertex). Therefore **every `best_count` reported by explore_families.py may be spuriously low**; only the `euclid_count` values and the exact checks in circles_exact.py are trustworthy. Also, the module runs its full exploration at import time (no `if __name__ == "__main__"` guard), so importing it hangs. Fix: reject candidates within a tolerance (say 1e-5) of any point of P, and guard the main block.

## 10. Suggested next steps for the computational team
1. Treat n = 8 as the first target: is 18 optimal? A Möbius-type enumeration (partial 3-linear spaces on 8 points with few blocks whose every derived structure is a realisable line arrangement) plus realisability checks would settle c(8); 56 triples, blocks of size 3–7.
2. For n = 5…7 the same approach is small enough for exhaustive SAT/enumeration; results would fill issue #297 (human-verified computation required for OEIS).
3. Re-derive Elliott's counting with Csima–Sawyer / exact small ordinary-line numbers to push the threshold 393 down — the most plausible "theory" contribution.

## Sources
- https://www.erdosproblems.com/506 ; https://www.erdosproblems.com/forum/thread/506 ; https://www.erdosproblems.com/tags/geometry
- https://raw.githubusercontent.com/google-deepmind/formal-conjectures/main/FormalConjectures/ErdosProblems/506.lean (copied to 506.lean)
- https://github.com/google-deepmind/formal-conjectures/issues/2128 ; https://github.com/google-deepmind/formal-conjectures/pull/4396 ; https://github.com/google-deepmind/formal-conjectures/issues/4952
- https://github.com/teorth/erdosproblems/issues/297 ; https://github.com/teorth/erdosproblems/wiki/AI-contributions-to-Erd%C5%91s-problems
- Elliott 1967: Acta Math. Acad. Sci. Hungar. 18, 181–188 (Springer landing page https://link.springer.com/doi/10.1007/BF01447430 surfaced by search)
- Purdy–Smith 2010: https://arxiv.org/abs/0907.0724 ; https://link.springer.com/article/10.1007/s00454-010-9270-3
- Bálintová–Bálint 1994: Acta Math. Hungar. 63, 283–289 (cited in the above and in arXiv:1607.06597, 1905.09639)
- Zhang 2011: https://link.springer.com/article/10.1007/s00454-010-9286-8
- Lin–Makhul–Mojarrad–Naslund–Swanepoel: https://arxiv.org/abs/1412.8314
- Lin–Makhul–Mojarrad–Schicho–Swanepoel–de Zeeuw 2018: https://arxiv.org/abs/1607.06597 ; https://link.springer.com/article/10.1007/s00454-017-9885-8
- Lin–Swanepoel 2019/2021: https://arxiv.org/abs/1905.09639
- Alon–Mellinger–Mubayi–Verstraëte: https://arxiv.org/abs/1007.4150
- Erdős 1961 problem list: https://renyi.hu/~p_erdos/1961-22.pdf (not fetched)
- OEIS (related, not this sequence): https://oeis.org/A003035 (orchard), https://oeis.org/A186704 (distinct distances)
