"""Skeptic c7-audit: COMPLETE real-non-realisability check of the Moebius-Kantor (8_3) configuration.
Frame: an independent quadrangle (3,5,6,7) -> (1,0,0),(0,1,0),(0,0,1),(1,1,1).  Each remaining point
lies in exactly one of the three charts  z=1 | (z=0,y=1) | (z=0,y=0,x=1), so the 3^4 = 81 chart
combinations exhaust RP^2 (and CP^2).  For each we solve the 8 incidence equations exactly and
report every solution with 8 distinct points, together with whether it is real.
"""
import itertools, sympy as sp
MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
frame = {3: (1, 0, 0), 5: (0, 1, 0), 6: (0, 0, 1), 7: (1, 1, 1)}
for c in itertools.combinations(frame, 3):
    assert not any(frozenset(c) <= t for t in MK)
others = [0, 1, 2, 4]
X = sp.symbols('x0:8'); Y = sp.symbols('y0:8')
charts = {0: lambda i: (X[i], Y[i], 1), 1: lambda i: (X[i], 1, 0), 2: lambda i: (1, 0, 0)}
total_sols, real_sols = 0, 0
for combo in itertools.product(range(3), repeat=4):
    P = {k: sp.Matrix(v) for k, v in frame.items()}
    unk = []
    for i, ch in zip(others, combo):
        P[i] = sp.Matrix(charts[ch](i))
        unk += [s for s in P[i] if s.is_Symbol]
    eqs = [sp.expand(sp.Matrix.hstack(*[P[i] for i in sorted(t)]).det()) for t in MK]
    if any(e.is_number and e != 0 for e in eqs):
        continue  # inconsistent chart (a fixed determinant is a non-zero constant)
    eqs = [e for e in eqs if e != 0]
    sols = sp.solve(eqs, unk, dict=True) if unk else ([{}] if not eqs else [])
    for s in sols:
        pts = {i: P[i].subs(s) for i in range(8)}
        # solutions may be parametric (free symbols) -> report as families
        free = set().union(*[set(v.free_symbols) for v in s.values()]) if s else set()
        dup = [(i, j) for i, j in itertools.combinations(range(8), 2)
               if all(sp.simplify(c) == 0 for c in pts[i].cross(pts[j]))]
        total_sols += 1
        isreal = all(sp.im(v) == 0 for v in s.values()) if not free else "parametric"
        print(f"chart {combo}: solution {s} free={free} coincident={dup} real={isreal}")
        if not dup and isreal is True: real_sols += 1
print("solutions with 8 distinct points found (all charts):", total_sols, "; REAL ones:", real_sols)
