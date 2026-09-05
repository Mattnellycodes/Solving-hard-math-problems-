# Solving Hard Math Problems

An attempt to make measurable, independently verifiable progress on an unsolved mathematics problem
using AI-orchestrated research: parallel hypothesis-testing agents, computational search, and adversarial
verification of every claim.

## The problem chosen: Erdős problem #506

*What is the minimum number c(n) of circles determined by n points in the plane, not all on one line
and not all on one circle?* (A circle is determined if it passes through at least three of the points;
collinear triples determine no circle.) Elliott (1967, corrected by Purdy–Smith) proved
c(n) = C(n−1,2) + 1 − ⌊(n−1)/2⌋ for n > 393. The problem was open for every 4 ≤ n ≤ 393, no small value
had ever been computed, and no OEIS sequence existed.

## Results (2026-09-05)

| n | c(n) | formula | status |
|---|------|---------|--------|
| 4 | 3 | 3 | elementary |
| 5 | 5 | 5 | proved (hand + enumeration) |
| 6 | **8** | 9 | **proved; formula fails** — triangle + feet of its altitudes |
| 7 | **11** | 13 | **proved; formula fails** — orthocentric system + its three diagonal points |
| 8 | **17** | 19 | **proved; formula fails** — beats Segre's classical cube example (18) |
| 9 | 25 | 25 | proved; independently re-verified (complete enumeration) |
| 10 | 33 | 33 | **proved; independently audited** (two enumerations reduce to one geometric structure, excluded by exact algebra) |
| 11–16 | ≤ f(n) | | no configuration below the formula found by two searches; conjectured equal |
| ≥ 11 | ? | | conjectured equal to the formula |

Sequence c(n), n = 4..10: **3, 5, 8, 11, 17, 25, 33** (all proved and independently verified). Full write-up: `experiments/506/RESULTS.md`; research note: `docs/NOTE_506.md`.

Certificates: exact integer coordinates for every upper bound (`experiments/506/circles_exact.py` recounts
them in rational arithmetic); computer-assisted lower bounds via a Möbius-plane reformulation
(circles = C(n,3) − D − ℓ), a CP-SAT / orderly enumeration of abstract block structures under
Sylvester–Gallai-derived constraints, and two short geometric lemmas. Each lower bound was re-derived
from scratch by the main agent and by two independent skeptic agents.

Explicit 8-point configuration with 17 circles: (0,0), (0,5), (0,10), (0,15), (3,6), (3,9), (5,5), (5,10).


## Second result: OEIS A067857 sign conjecture refuted (2026-09-05)

The problem-selection survey ranked first a conjecture of Robert Israel (OEIS A067857, 2015; formalised as
research-open in google-deepmind/formal-conjectures): with a(n) defined by Σ_{k|n} a(k)/k! = H_n, the claim was
that a(n) < 0 **iff** ω(n) is odd and ≥ 3. It is false in both directions:

- a(479#) < 0 although ω(479#) = 92 is even; a(487#) > 0 although ω(487#) = 93 is odd (p₉₂ = 479, p₉₃ = 487);
- **Theorem:** sign a(p_k#) = (−1)^{k+1} for every k ≥ 92 (uniform monotone bound, one finite exact-rational
  computation UB(92) < 0; no Mertens-type input);
- corrected sign law: for squarefree n with ω(n) ≥ 2, a(n)/n! = (−1)^{ω(n)} G(n) with G an explicit alternating
  sum of Euler products plus a bounded error, so the sign flips once ∏_{p|n}(1+1/p) exceeds ≈ 6.74;
- the conjecture holds for all n ≤ 4000 and every counterexample has ω(n) ≥ 52; minimality of 479# is open.

Certificates: five independent implementations (mpmath interval arithmetic, exact rationals, a pure-Python
fixed-point checker, MPFR directed rounding, and the main agent's own check) agree; both hostile audits
returned CONFIRMED. See `experiments/A067857/RESULTS.md`, `proof/PROOF.md`, `certificate/`.

## Layout

- `docs/JOURNAL.md` — chronological research journal (every experiment, every result, including negative ones)
- `docs/PROBLEM_506.md` — problem dossier and hypothesis table
- `docs/ENVIRONMENT.md` — compute, network, and data sources
- `experiments/506/theory/` — theory agent's report and enumeration code
- `experiments/506/literature/` — literature report
- `experiments/506/verify_independent/` — main agent's independent verification (CP-SAT model, certificates)
- `experiments/506/skeptics/` — adversarial verification by independent agents
- `experiments/A067857/` — refutation of the A067857 sign conjecture (proof, certificates, audits)
- `scripts/` — shared tooling
