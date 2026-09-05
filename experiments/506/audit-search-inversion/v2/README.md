# audit-search-inversion / v2 — independent re-derivation (second audit pass)

All code in this directory was written from scratch for this audit pass; nothing from
`search-inversion/` (or from the earlier audit files in the parent directory) is used as evidence.

| file | what it does | result |
|---|---|---|
| `count506.py` | exact planar counter: canonical integer key (a,b,c,d) of the generalised circle a(x²+y²)+bx+cy+d=0 through each triple (3×3 minors), integer-only; recomputes every block's membership independently and checks the triple partition | — |
| `records.py` | recounts the 8 antipodal records n=9..16 from the coordinates quoted in the report | all = f(n), non-degenerate, structure = one (n−1)-circle + ⌊(n−1)/2⌋ 3-point lines + 3-point circles |
| `sanity.py` | validates `count506` on the certified n=6,7,8 sets of `theory/REPORT.md` (8, 11, 17) and against `circles_exact.py` on 300 random sets | all agree |
| `s2lat_check.py` | the two lattice-sphere n=9 sets: exact rational inversion of R³ about the sphere point Q, then a 3-D coplanar circle/line count (circumcentre + r² keys) — a different computation from plane counting on the sphere | S2lat-5 and S2lat-6: 25 circles, 4 lines, non-degenerate, for the certified centres (2,0,−1) / (−2,−1,1) and for the "description" centres (−2,−1,0) / (−2,−1,−1) |
| `budget.py` | recomputes the D+l thresholds C(n−1,3)+⌊(n−1)/2⌋ and the all-4-block caps 3n·t3(n−1)/4 + t3(n) | 60,88,125,170,226,292,371,462 and 57,87,115,163 as claimed |
| `pentagon.py` | two concentric regular pentagons at 4 radius ratios × 2 orientations (mpmath, 50 digits) | D = 78 in every case (2 five-circles + 20 four-blocks); pentagram ratio gives 5 four-point lines and 37 circles |
| `sphsearch.py` | exhaustive subset + projection-centre search on x²+y²+z²=N, every centre with ≥2 blocks enumerated exactly (rational points or irrational conjugate pairs); rational-centre hits re-verified by the projection counter | see table below |

Orchard numbers t3(3..12) = 1,1,2,4,6,7,10,12,16,19 checked against OEIS A003035 (Burr–Grünbaum–Sloane 1974).

## Lattice-sphere table (exhaustive; run logs run*.log)
(filled in below)
Runs relaunched at ~20:58 with the corrected survivor bound (L <= floor(C(n,2)/3)); at the time this
audit pass was forced to report, the N=5 and N=6 exhaustive runs (n=9..16) and the chained N=10, N=11
runs had not yet printed results (logs run5.log, run6.log, run10.log, run11.log). The earlier audit pass
in the parent directory (run_A.log, run_B.log, run_N5.log, run_N6.log) reports, with its own code:
N=5: n=9 25, n=10 39, n=11 45, n=12 61; N=6: n=9 25, n=10 45, n=11 59, n=12 78; N=10: n=9 25, n=10 45;
N=11: n=9 25, n=10 39 -- all agreeing with the search-inversion table where it is OPTIMAL/claimed.
Also verified here (verify_results2.py): every S2lat row of results2.jsonl reaching f(n) for n=10..13
(17 rows on N=14,17,26,29,38) recounts exactly to f(n) and is of antipodal Moebius type.
