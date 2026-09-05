# Erdős problem #506 — minimum number of circles determined by n points

**Statement.** Let P be a set of n points in the plane, not all on one circle and not all on one line.
A circle is *determined* by P if it passes through at least three points of P. What is the minimum
number c(n) of circles determined by such a set?

**Known.** Elliott (1967), corrected by Purdy–Smith, proved for n > 393 that
c(n) = C(n−1,2) + 1 − ⌊(n−1)/2⌋, attained by n−1 points on a circle in antipodal pairs plus the
centre (the antipodal triples through the centre are collinear and determine no circle).
The problem is open for 4 ≤ n ≤ 393. Segre observed that for n = 8 a projected cube beats the
uncorrected bound C(7,2) = 21.

**Möbius reformulation (used by all code here).** Let B(P) be the set of blocks: circles *or lines*
through ≥ 3 points of P. Inverting about a point O ∉ P maps blocks through O to lines, so the
Euclidean circle count of the inverted set is |B(P)| − deg(O), where deg(O) is the number of blocks
through O. Taking O = ∞ recovers circles(P) = |B(P)| − #lines. Hence the best count reachable from
the Möbius configuration P is |B(P)| − max_O deg(O).

**First result (2026-09-05, verified exactly).** Two concentric axis-aligned squares
{(±1,0),(0,±1),(±2,0),(0,±2)} determine exactly 18 circles (8 collinear triples on the two axes,
2 circumcircles, 8 further 4-point circles, 8 ordinary circles), one fewer than the large-n
formula value 19. Whether 18 is optimal for n = 8, and the values for n = 9 … 393, are open.

Files:
- `circles_exact.py` — exact (rational) circle counter; certificate checker.
- `explore_families.py` — numeric exploration of structured families with inversion-centre optimisation.
