"""Sanity check of the Miquel/bundle filters and of the skeleton-3 geometry: two concentric regular
pentagons with parallel sides (R = -lambda * S) have exactly the 20 concyclic 'parallel-chord' quadruples."""
import math, itertools, sys
import numpy as np
sys.path.insert(0, '.')
from realise import structure_of
from filters import miquel_violations, bundle_violations, full_report
for lam in (0.5, math.cos(math.pi/5), 0.7):
    S = [(math.cos(math.pi/2 + 2*math.pi*k/5), math.sin(math.pi/2 + 2*math.pi*k/5)) for k in range(5)]
    R = [(-lam*x, -lam*y) for x, y in S]
    P = np.array(S + R)
    F, L, co = structure_of(P, tol=1e-9)
    print(f"lambda={lam:.4f}: rich blocks {sorted(len(B) for B in F)}, lines {[sorted(l) for l in L]}")
    print("   Miquel violations:", len(miquel_violations(F, L)), " bundle violations:", len(bundle_violations(F, L)),
          " count = 120 -", sum(math.comb(len(B),3)-1 for B in F), "-", len(L), "=", 120 - sum(math.comb(len(B),3)-1 for B in F) - len(L))
    F4 = [sorted(B) for B in F if len(B) == 4]
    print("   4-blocks:", F4)
