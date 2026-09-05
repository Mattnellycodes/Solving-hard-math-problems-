# Independent verification log (main agent, 2026-09-05)

All code in this directory was written from scratch by the main session agent, without reusing the
theory agent's enumeration code, to check the claims in `../theory/REPORT.md`.

## Constructions (exact rational arithmetic, `../circles_exact.py`)
| n | points | circles | formula f(n) |
|---|--------|---------|--------------|
| 6 | (0,0),(20,0),(5,15),(5,0),(2,6),(10,10) | **8** | 9 |
| 7 | the six above + (5,5) | **11** | 13 |
| 8 | (0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10) | **17** | 19 |
| 8 | two concentric squares (±1,0),(0,±1),(±2,0),(0,±2) | 18 | 19 |
All non-degenerate (not all concyclic / collinear).

## Combinatorial relaxation (`cpsat_relaxation.py`, OR-tools CP-SAT)
Variables: rich blocks (subsets of size 4..n−1) and lines (subsets of size 3..n−1). Constraints: two
rich blocks share ≤ 2 points; a line of size ≥ 4 is a rich block; a 3-line is an uncovered triple; two
lines share ≤ 1 point; derived Sylvester–Gallai cap at every point; line cap. Objective max D + ℓ.
Run with the weakest safe caps (`sg`: o(m) = 1) and with the known table (`table`).

| n | opt D+ℓ (sg) | opt D+ℓ (table) | ⇒ circles ≥ | structures at the optimum |
|---|---|---|---|---|
| 6 | 13 | 13 | 7 | all: 3 four-blocks = complements of a perfect matching + 4 transversal 3-lines (complete quadrilateral); 30 labelled copies |
| 7 | 27 | 27 | 8 | all structures with count ≤ 10: the 7 four-blocks = Fano complements, with 4, 5 or 6 lines |
| 8 | 39 | 39 | **17** | — (so c(8) ≥ 17 needs only Sylvester–Gallai) |
| 9 | 65 | (running) | 19 | optimum has 18 four-blocks, every point in 8 → derived structure is (8_3) |

Hence **c(8) = 17 is fully verified**: lower bound from the relaxation with SG-only caps, upper bound from
the integer construction.

## Geometric steps
* **Angle Lemma** (`angle_lemma_poly.py`): with L1 the x-axis and L_i: y = m_i x + c_i, the three
  concyclicity determinants factor as (degeneracy factors: P_ij = P_ik or three lines concurrent) ×
  A_i with A1 = m2m3m4 + m2 + m3 − m4, A2 = m2m3m4 + m2 − m3 + m4, A3 = m2m3m4 − m2 + m3 + m4, and
  A1 − A2 = 2(m3 − m4) ≠ 0. So the three diagonal quadruples of a complete quadrilateral are never all
  concyclic. (Also checked the directed-angle proof by hand.)
* **n = 6**: in every relaxation-optimal structure the 4 lines form a complete quadrilateral whose
  opposite-vertex pairs are exactly the matching pairs, so the three 4-blocks are the diagonal
  quadruples ⇒ impossible ⇒ **c(6) = 8**.
* **n = 7**: in the Fano-complement structure, at any point p the four blocks through p pairwise share
  exactly one further point, so after inversion about p they are four pairwise non-parallel,
  non-concurrent lines (a complete quadrilateral on the other 6 points) and the three blocks avoiding p
  are its diagonal quadruples ⇒ impossible ⇒ **c(7) = 11**.
* **(8_3)** (`mk_unrealisable.py`, `mk_unrealisable2.py`): 840 families of 8 pairwise ≤1-intersecting
  triples on 8 points = one S_8-orbit of size 8!/48 (Aut of Möbius–Kantor has order 48), so every such
  family is Möbius–Kantor. With frame p0,p1,p2,p5 and forced coordinates, realisability requires
  s = 1 and t² − t + 1 = 0 — no real solution. (My first attempt used a wrong chart and is discarded.)
* **n = 9** (`cpsat_n9.py`, `cpsat_n9_sym.py`): adds "≤ 7 four-blocks through any point" (valid by (8_3))
  and applies an independent hereditary-SG check to every structure with ≤ 24 circles. The symmetry-free
  run hit its 2-hour limit (status FEASIBLE, incomplete). The symmetry-broken run — justified because
  D + ℓ ≥ 60 with b₄ ≤ 15 and ℓ ≤ 11 forces a block of size ≥ 5, which we relabel to contain {0,1,2,3,4} —
  finished in 197 s with status OPTIMAL (complete): 720 labelled solutions, all one isomorphism class
  (two 5-blocks sharing a point + twelve 4-blocks + 6 lines, count 24), and **all 720 fail hereditary SG**
  (a derived Fano plane at a degree-7 point). Hence no structure with ≤ 24 circles is realisable, and with
  the antipodal configuration (25) we get **c(9) = 25**, independently verified with SG-only caps.
  Dependencies: Sylvester–Gallai (for the caps and the hereditary check) and the (8_3) non-realisability
  lemma (`mk_unrealisable2.py`). No orchard or Kelly–Moser values are used.
