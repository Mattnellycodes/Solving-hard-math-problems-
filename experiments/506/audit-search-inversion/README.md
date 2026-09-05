# audit-search-inversion — independent re-derivation tools

All code here was written from scratch for the audit (no reuse of search-inversion/ code).

* `audit_count.py`   exact planar counter (Fractions, circumcentre keys); recounts the 8 antipodal records.
* `det_count.py`     second planar recount via the 4x4 concyclicity determinant (integer arithmetic only).
* `sphere_audit.py`  exact 3D plane counting on integer lattice spheres; audits the two S2lat n=9 sets.
* `sphere_search.py` exhaustive subset + projection-centre search on x^2+y^2+z^2=N
                     (`python3 sphere_search.py N n:target [n:target ...]`); logs in run_N5.log, run_N6.log, run_A.log, run_B.log.
* `project_check.py` cross-check of search hits by numerical centre + stereographic projection + 40-digit planar count.
* `pent_check.py`    block structure of two concentric regular pentagons (D = 78).
