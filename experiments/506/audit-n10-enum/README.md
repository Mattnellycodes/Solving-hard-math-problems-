# audit-n10-enum: independent re-derivation of the n10-enum result (own code only)

Scripts (all written from scratch here; nothing imported from other agents' directories):
* `mob.py`      – model, brute-force S_n permutation table, stage 1 (blocks >= 5, canonical = lex-min sorted
                  bitmask list under all of S_n), stage 2 (4-blocks, canonical under Aut(F_big)), exact line-set
                  enumeration, brute-force canonical forms of full structures.
* `validate.py` – V1: stage1+stage2 reproduce the brute-force isomorphism classes of all labelled rich-block
                  families for n = 6, 7 (each class once); V2: stage-2 canonicity test vs brute force (n = 7, 8).
* `run.py n target mode [--all4]` – full pipeline (big-block case; `--all4` adds the empty skeleton with G = S_n).
* `all4.py`     – n = 10 all-4-block case via a maximum-degree point and its derived PSTS(9, d0), d0 in {10, 11}.
* `cpsat_b4_25.py`, `cpsat_all4.py` – independent OR-tools CP-SAT checks of the all-4-block branches.
* `checks.py`   – tiers A (hereditary SG), B ((8_3)), K (Kelly–Moser), C (orchard / o-tables) on the 11 derived
                  real point sets of Q = P + {inf}, for every (F, L).
* `crosscheck.py` – brute-force S_10 canonical forms: my classes vs the other agent's JSON.
* `aut_and_mk.py` – |Aut(F)| of the survivors; uniqueness of the 8-triple (8_3) systems.

Logs: `validate.log`, `n8_t18.log`, `n9_t25.log`, `n10_t32_sg.log`, `n10_t33_sg.log`, `all4_n10_t32.log`,
`cpsat_all4.log`, `checks_n10_t32.log`, `checks_all4.log`, `checks_n9.log`, `crosscheck_n10.log`,
`crosscheck_all4.log`, `aut_and_mk.log`; JSON outputs `out_*.json`.
