# audit-n10-cpsat — independent audit of the n10-cpsat (v2) result (Erdős #506, n = 10)

All code here was written from scratch by the auditor; nothing is imported from n10-cpsat, n10-continue,
theory or verify_independent (their outputs are only compared against).

Files: `common.py` (structures, isomorphism classes, real-linear-space filters: SG cap, hereditary SG,
(8_3), orchard caps, Miquel closure), `miquel_check.py` (cross-ratio identity), `skeletons.py`
(families of blocks of size >= 5 -> `skeletons_audit.json`, 22 classes), `enum_f4.py` (stage 1, CP-SAT
all-solution enumeration of 4-block families per skeleton with own lex-leader symmetry breaking;
`--nolex`, `--nomk`, `--t3cap`, `--first`, `--seed=`, index -1 = no big block), `lines.py` (stage 2,
exhaustive line sets), `post.py` (filters), `two5_exact.py` + `snf.py` (exact Möbius normal-form /
Smith-normal-form / Gröbner analysis for two disjoint 5-blocks), `sanity_antipodal.py`,
`sanity_two5.py`, `mk_check.py` ((8_3) lemma), `match_v2.py` (locates the v2 examples), `summarize.py`.
Outputs: `f4_*.json/.log`, `fl_*.json/.log`, `post_*.log`, `surv_*.json`, `two5_*.log`, `*.log`.
