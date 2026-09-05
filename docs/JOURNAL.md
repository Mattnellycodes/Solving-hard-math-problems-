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
