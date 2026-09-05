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
