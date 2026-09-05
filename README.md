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
| 9 | 25 | 25 | proved by the theory agent; independent re-verification in progress |
| ≥ 10 | ? | | conjectured equal to the formula |

Certificates: exact integer coordinates for every upper bound (`experiments/506/circles_exact.py` recounts
them in rational arithmetic); computer-assisted lower bounds via a Möbius-plane reformulation
(circles = C(n,3) − D − ℓ), a CP-SAT / orderly enumeration of abstract block structures under
Sylvester–Gallai-derived constraints, and two short geometric lemmas. Each lower bound was re-derived
from scratch by the main agent and by two independent skeptic agents.

Explicit 8-point configuration with 17 circles: (0,0), (0,5), (0,10), (0,15), (3,6), (3,9), (5,5), (5,10).

## Layout

- `docs/JOURNAL.md` — chronological research journal (every experiment, every result, including negative ones)
- `docs/PROBLEM_506.md` — problem dossier and hypothesis table
- `docs/ENVIRONMENT.md` — compute, network, and data sources
- `experiments/506/theory/` — theory agent's report and enumeration code
- `experiments/506/literature/` — literature report
- `experiments/506/verify_independent/` — main agent's independent verification (CP-SAT model, certificates)
- `experiments/506/skeptics/` — adversarial verification by independent agents
- `scripts/` — shared tooling
