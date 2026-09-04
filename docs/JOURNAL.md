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

