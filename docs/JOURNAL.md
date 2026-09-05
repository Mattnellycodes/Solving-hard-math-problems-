# Research Journal

All times UTC. Every entry records what was tried, what happened, and what it means.
Negative results are recorded with the same care as positive ones.

---

## 2026-09-04 — Session start

**Goal.** Pick an unsolved problem where computation plus careful reasoning has a real chance
of producing a verifiable improvement (a new bound, a new extremal construction, or a settled
small case), then attack it with parallel hypothesis-testing workflows.

**Environment.** 4 CPUs, 15 GB RAM, Python 3.11 with sympy/numpy/scipy/networkx/z3/pysat/ortools/gmpy2/numba.

**Selection criteria.**
1. Precisely stated, finite or computationally checkable objective (no "prove X for all n" as the only goal).
2. Known best bounds are documented and recent enough to be trusted.
3. Improvement is verifiable by an independent program (a certificate).
4. Not a famous problem with decades of expert effort (Riemann, Goldbach, Collatz, twin primes are out).
5. Ideally from a curated list (Erdős problems, OEIS, Guy's UPINT) so "open" status is authoritative.

**Reconnaissance (see `docs/ENVIRONMENT.md`).** Most math reference sites are blocked by the egress
proxy; GitHub and web search work. Cloned the Erdős-problems community database, the formal-conjectures
Lean repository, and the LLM-attack archive locally. The database exposes three especially interesting
status classes: 9 problems marked *decidable* (open but reduced to a finite computation), 25 *falsifiable*,
7 *verifiable*.

The nine decidable problems: 19 (chromatic number, $500), 475 (distinct partial sums of subsets of F_p),
506 (minimum number of circles determined by n points, open for n ≤ 393), 547/551/556 (Ramsey theory),
580 (Loebl–Komlós–Sós n/2 conjecture, proved for large n by Zhao), 742 (Murty–Simon diameter-2-critical
conjecture, verified n ≤ 24 and n = 26 by Fan, all large n by Füredi), 848 (ab+1 never squarefree;
asymptotic case proved by Sawhney).

**Workflow 1 launched — problem selection** (`select-unsolved-problem`, run wf_ce7e80cd-c21).
Ten scouts with distinct lenses (decidable / falsifiable ×2 / verifiable / open-with-OEIS / open-with-explicit-bounds /
OEIS conjectures / other formal-conjectures collections / classical computational records / AI-search frontier),
then two referee lenses per shortlisted candidate (already-done check; tractability check with small calibration
experiments allowed), then a synthesis agent producing a ranked list and a survey report.


**Push blocked.** `git push` and the GitHub API both return 403 for this repository ("Claude doesn't have
GitHub access ... / Resource not accessible by integration"). The Claude GitHub App appears not to be
installed for the `Mattnellycodes` account or not granted this repo. Committing locally; will retry pushes.

**Pre-research on a-priori favourite, Erdős #475 (Graham's rearrangement conjecture).**
Status has moved fast in 2024–2026:
- Costa–Della Fiore–Ollis–Rovner-Frydman (EJC 2022): true for |A| ≤ 12, all primes, via the polynomial method
  (Alon's Combinatorial Nullstellensatz; polynomial of degree k²−k−1 in k variables).
- Kravitz (Integers 2024): |A| ≤ log p / log log p. Bedert–Kravitz (Israel J. Math 2025): |A| ≤ exp((log p)^{1/4}).
- Costa–Della Fiore (arXiv 2602.19989, Feb 2026): |A| ≤ exp(c (log p)^{1/3}).
- **Pham–Sauermann (arXiv 2602.15797, Feb 2026): |A| ≤ p^{1−α} for |A| large given α; with earlier results this
  resolves the conjecture for all sufficiently large primes.** So what remains is small primes (ineffective
  threshold) and the non-prime cyclic / general abelian case ("Graham conjecture on small sets in abelian groups",
  arXiv 2603.20961, Mar 2026). Extending the all-p verification from |A| ≤ 12 to 13 is well-defined but of
  reduced interest now. Left to the judges.

**Pre-research on Erdős #742 (Murty–Simon).** Fan's n ≤ 24 is a proof, not an enumeration. Kirchweger–Szeider
(2024, SAT modulo symmetries) enumerated all diameter-2-critical graphs up to 13 vertices. Extending either
route (proof for n = 25, or enumeration to n = 14) is far beyond 4 CPUs. Deprioritized.

**Pre-research on Erdős #506 (minimum number of circles).** Sanity-checked the large-n formula
C(n−1,2) + 1 − ⌊(n−1)/2⌋: it is attained by n−1 points on a circle in antipodal pairs plus the centre (the
⌊(n−1)/2⌋ collinear triples through the centre determine no circle). For n = 4 the minimum is 3 (matches).
Open question for n ≤ 393 is whether some configuration with many concyclic quadruples beats this. This is a
"find a construction" problem — attractive because any construction is a certificate.

## 2026-09-05

**Survey workflow slowed by the 2-agent concurrency cap** (two scouts finished in ~40 min). Stopped it,
cut the judge phase to one merged referee for the top 8 candidates, and resumed (cached scouts reused).
Scout findings so far: Erdős #475 (Graham's rearrangement conjecture) is now verified for |A| ≤ 20
(Costa–Della Fiore–Fontana–Vena, arXiv 2603.20961) and settled for all large p (Pham–Sauermann), so the
"extend |A| ≤ 12 to 13" idea is dead. #23 (Ferudun, arXiv 2606.28041, June 2026) proved a(5n) = n² for n ≤ 40
by flag algebras; extending is beyond 4 CPUs.

**Independent lead: Erdős #506 (minimum number of circles determined by n points).** While waiting I worked
this problem by hand and by computer (details in `experiments/506/README.md`):
- Möbius reformulation: circles(P) = |B(P)| − deg(O) where B(P) are the circles-or-lines through ≥ 3 points and
  O is the inversion centre (O = ∞ gives the Euclidean count, deg(∞) = number of ≥3-point lines). Equivalently
  circles(P) = C(n,3) − Σ_{circle blocks}(C(k,3) − 1) − Σ_{line blocks} C(k,3).
- **n = 8: two concentric axis-aligned squares determine exactly 18 circles** (verified with exact rational
  arithmetic, `experiments/506/circles_exact.py`), one below the large-n formula value 19. This is probably
  Segre's "projected cube" (perspective projection along a 4-fold axis) — to be confirmed by the literature agent.
- Hand analysis: any configuration with two or more points off the main circle costs Θ(n²) extra circles, so
  exceptions to the formula can only live at small n (block-size-4 designs beat n²/2 only for n ≲ 15).
  SQS(8) is not realisable (its derived design is the Fano plane); the 7-block "anti-Fano" quadruple system on
  7 points would tie, not beat, the formula at n = 7 (two of its blocks share two points, so at most one can be a line).
- Not previously tabulated as far as web search shows: Bálintová–Bálint (1994), Zhang (2011) and
  Lin–Makhul–Mojarrad–Naslund–Swanepoel (2017) concern *ordinary* circles; the erdosproblems database lists
  no OEIS sequence for #506 ("possible").

**Workflow 2 launched — Erdős #506 deep dive** (`erdos-506-deep-dive`, run wf_cdf13f56-0d0): literature scout,
lower-bound theorist (rigorous c(n) for n = 5,6,7, attempt 8, computer-assisted plan), four construction
searchers (symmetric orbits with continuous parameters; subset/local search over rich universes incl. exhaustive
8-subsets of 4×4 and 5×5 grids; classical configurations and stereographic projections of polyhedra; hybrid
near-pencil families), an independent exact verifier per searcher, and a synthesis agent.

## 2026-09-05 (continued) — verification of the #506 results

The usage limit interrupted both workflows at ~01:50 UTC (all 10 scouts and the #506 literature and
theory agents had finished; referees, searchers and synthesis had not). Resumed at 05:50 UTC.

**Theory agent's claims** (`experiments/506/theory/REPORT.md`): c(5)=5, **c(6)=8, c(7)=11, c(8)=17**
(all three below the large-n formula), c(9)=25 (= formula); explicit integer constructions; lower bounds
by enumerating abstract Möbius block structures under Sylvester–Gallai-derived necessary conditions,
finished by an "Angle Lemma" (n=6,7) and by non-realisability of the Fano and Möbius–Kantor configurations
(n=9). Conjecture: c(n) = f(n) for all n ≥ 9.

**Literature agent** (`experiments/506/literature/REPORT.md`): no exact value of c(n) for any 4 ≤ n ≤ 393
exists in the literature; no OEIS sequence (teorth/erdosproblems issue #297 asks for exactly this
computation). Segre's cube example is the two-concentric-squares configuration (18 circles), so my
18 was a rediscovery; 17 is new.

**Independent verification by the main agent** (`experiments/506/verify_independent/`):
- Constructions recounted exactly: 8, 11, 17 circles. ✔
- Own CP-SAT model of the relaxation: optimum D+ℓ = 13, 27, 39 for n = 6, 7, 8 (identical with SG-only
  and with the Kelly–Moser table caps). So **c(8) ≥ 17 follows from Sylvester–Gallai alone; c(8) = 17 is
  established.** For n = 6, 7 the relaxation leaves exactly the Angle-Lemma structures (verified by
  enumeration: matching-complement blocks + complete quadrilateral; Fano-complement blocks).
- Angle Lemma: my own factorisation of the three concyclicity determinants gives essential factors with
  A1 − A2 = 2(m3 − m4); hand-checked the directed-angle proof too. ✔ ⇒ **c(6) = 8, c(7) = 11 established.**
- (8_3): all 840 admissible 8-triple systems on 8 points form one S_8-orbit (Möbius–Kantor, |Aut| = 48);
  realisation forces t² − t + 1 = 0 ⇒ not realisable over ℝ. ✔ (First attempt used a wrong chart; fixed.)
- n = 9: my SG-only relaxation gives only ≥ 19 (optimum uses 18 four-blocks with every point in 8, i.e.
  (8_3) derived structures). Added the (8_3)-derived constraint and an independent hereditary-SG checker;
  enumeration of all structures with ≤ 24 circles is running.

**Adversarial workflow** (`erdos-506-adversarial-verify`, run wf_c96534f2-39f): two skeptics per claim
(construction attack / proof audit with own code) plus adjudication — running. First verdict in: c(6)=8
CONFIRMED by an independent enumeration.

**GitHub push now works** (the app was evidently granted access); branch pushed.

## 2026-09-05 10:45 UTC — skeptic verdicts and n = 9 status

**Adversarial verification (6 of 10 skeptic reports in; the rest were cut by the usage limit and are re-running):**
| claim | construction-attack | proof-audit |
|---|---|---|
| c(6) = 8 | CONFIRMED — own enumeration with only the intersection axioms gives the unique 7-circle candidate; own slope algebra kills it; exact brute force over ~10⁹ six-subsets of a 7×7 grid and a 97-point rational set: minimum 8 | CONFIRMED — own brute-force + z3 enumeration, own two-chart Angle-Lemma certificate |
| c(7) = 11 | CONFIRMED — every ≤10-circle structure is the Fano line system (impossible) or the (7,4,2) biplane, shown unrealisable by an all-chart Gröbner analysis | CONFIRMED — own labelled enumeration + own CP-SAT; Angle Lemma re-proved by hand and by Gröbner |
| c(8) = 17 | CONFIRMED — CP-SAT optimum D+ℓ = 39 and a complete DFS (393,328 nodes) with the cube structure as the unique ≤17 class; searches over grids/polyhedra found nothing below 17 | CONFIRMED — own CP-SAT (different formulation, both cap regimes) and own orderly enumeration |

Useful remarks from the skeptics: the algebraic Angle-Lemma certificates that fix L1 horizontal with finite
slopes miss the chart where another line is perpendicular to L1; the synthetic directed-angle proof (which I
checked by hand) is complete, and the c6-audit skeptic supplied a two-chart certificate. Also: the problem's
non-degeneracy convention ("not all on a line" vs "no three collinear") is ambiguous on erdosproblems.com; all
our results use the Lean/Elliott convention (points not all on one line or one circle; collinear triples
determine no circle), and the constructions for n = 6, 7, 8 do contain collinear triples.

**n = 9.** My symmetry-free enumeration (SG-only caps + "≤ 7 four-blocks through a point", the latter justified by the
Möbius–Kantor lemma) found, before its 2-hour limit, only the theory agent's candidate (two 5-blocks through a
common point + twelve 4-blocks + 6 lines, count 24), and my hereditary-SG checker kills it (a derived Fano plane at
a degree-7 point). The run was cut off (status FEASIBLE), so completeness is not yet certified; re-running with
symmetry breaking (some block of size ≥ 5 must exist, so fix one to contain {0,1,2,3,4}).

**Survey.** All eight referees finished; the synthesis agent is re-running. The referees' top survivor is a
refutation of Robert Israel's sign conjecture for OEIS A067857 at primorials with ≥ 92 prime factors (with a
uniform bound argument) — a second promising target if time permits.

**11:00 UTC — c(9) = 25 independently confirmed.** Symmetry-broken complete enumeration
(`experiments/506/verify_independent/cpsat_n9_sym.py`, status OPTIMAL, 197 s): the only abstract structure with
≤ 24 circles compatible with Sylvester–Gallai caps and the (8_3) lemma is the two-5-block structure, and all
720 labelled copies fail hereditary Sylvester–Gallai. Together with the antipodal configuration this gives
c(9) = 25 = f(9), using only Sylvester–Gallai and the Möbius–Kantor lemma. So the exceptions to the
Elliott–Purdy–Smith formula among n ≤ 9 are exactly n = 6, 7, 8.

**11:20 UTC — survey synthesis complete** (`docs/PROBLEM_SELECTION.md`). 91 candidates, 8 refereed. Ranking:
1. OEIS A067857 sign conjecture (Robert Israel 2015; formalised as research-open in formal-conjectures):
   three independent scout/referee computations find it FALSE at the primorials 479# (ω = 92) and 487# (ω = 93),
   and a uniform monotone bound shows failure for every primorial with ≥ 92 prime factors. Cheap, fully
   certifiable — adopted as a second target.
2. Exact values of A341822(n) (Gowers–Long 2-increasing sequences of triples), n = 8.
3–5. Spencer constant C10c, ring-loading constant (AlphaEvolve P61), covering-design records.
Amusingly, Erdős #506 ranked 17th in the scouts' scoring ("not advanced" to referees) — a reminder that
tractability is hard to judge before trying.
