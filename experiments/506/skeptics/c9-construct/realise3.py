"""Numerical realisation attack, v3 (v2's residuals were not scale-invariant and accepted collapsed or
escaping points).  Residuals now: for a Mobius block, every 4-subset (z1..z4) as complex numbers must have
a real cross-ratio cr=(z1-z3)(z2-z4)/((z1-z4)(z2-z3)); residual = Im(cr)/(1+|cr|) -- invariant under
similarity and scale.  For a line structure, every triple gives residual = sin(angle) = det/(|pq||pr|).
Accepted realisation: residual < 1e-10, all pairwise distances >= SEP, all |coords| <= BOUND, and an
independent structural recount (tolerance clustering of the actual blocks) confirms every required block.
Complex-coordinate runs (cross-ratio residual = |Im| replaced by the full complex equation cr - conj(cr)
does not make sense over C, so for complex runs we use the determinant equations normalised by the row
norms) serve as positive controls for Mobius-Kantor, which IS realisable over C."""
import numpy as np, itertools, sys, time
from scipy.optimize import least_squares
from closure_search import block_vectors, TRIPLES
rng = np.random.default_rng(7)
SEP, BOUND = 0.05, 30.0

def unpack(X, n, fixed, cplx):
    P = np.zeros((n, 2), dtype=complex if cplx else float)
    for i, v in fixed.items(): P[i] = v
    free = [i for i in range(n) if i not in fixed]
    if cplx: X = X[:2 * len(free)] + 1j * X[2 * len(free):]
    P[free] = X.reshape(-1, 2)
    return P

def mob_res(X, blocks, fixed, n, cplx):
    P = unpack(X, n, fixed, cplx); z = P[:, 0] + 1j * P[:, 1] if not cplx else None
    out = []
    for b in blocks:
        for q in itertools.combinations(sorted(b), 4):
            if not cplx:
                z1, z2, z3, z4 = z[list(q)]
                cr = (z1 - z3) * (z2 - z4) / ((z1 - z4) * (z2 - z3) + 1e-300)
                out.append(cr.imag / (1 + abs(cr)))
            else:
                rows = np.array([[P[i, 0] ** 2 + P[i, 1] ** 2, P[i, 0], P[i, 1], 1.0] for i in q])
                nrm = np.prod(np.linalg.norm(rows, axis=1))
                out.append(np.linalg.det(rows) / nrm)
    out = np.array(out)
    return np.concatenate([out.real, out.imag]) if cplx else out

def line_res(X, lines, fixed, m, cplx):
    P = unpack(X, m, fixed, cplx); out = []
    for l in lines:
        for (i, j, k) in itertools.combinations(sorted(l), 3):
            u = P[j] - P[i]; v = P[k] - P[i]
            det = u[0] * v[1] - u[1] * v[0]
            nu = np.sqrt(u[0] ** 2 + u[1] ** 2 + 0j); nv = np.sqrt(v[0] ** 2 + v[1] ** 2 + 0j)
            out.append(det / (nu * nv + 1e-300))
    out = np.array(out)
    return np.concatenate([out.real, out.imag]) if cplx else out.real

def structure_ok(P, blocks, kind):
    """independent check: recompute the blocks of the real point set P by tolerance clustering and test
    that every required block is contained in a computed block."""
    P = np.array(P, float); n = len(P)
    if kind == "mob":
        found = {}
        for (i, j, k) in itertools.combinations(range(n), 3):
            rows = np.array([[P[t, 0] ** 2 + P[t, 1] ** 2, P[t, 0], P[t, 1], 1.0] for t in (i, j, k)])
            v = np.array([np.linalg.det(rows[:, [1, 2, 3]]), -np.linalg.det(rows[:, [0, 2, 3]]), np.linalg.det(rows[:, [0, 1, 3]]), -np.linalg.det(rows[:, [0, 1, 2]])])
            v = v / np.linalg.norm(v); v = v * np.sign(v[np.argmax(np.abs(v) > 1e-9)])
            found.setdefault(tuple(np.round(v / 1e-7).astype(np.int64)), set()).update((i, j, k))
    else:
        found = {}
        for (i, j) in itertools.combinations(range(n), 2):
            d = P[j] - P[i]; nrm = np.array([-d[1], d[0]]) / np.linalg.norm(d); c = nrm @ P[i]
            if nrm[0] < 0 or (abs(nrm[0]) < 1e-9 and nrm[1] < 0): nrm, c = -nrm, -c
            found.setdefault(tuple(np.round(np.array([nrm[0], nrm[1], c]) / 1e-7).astype(np.int64)), set()).update((i, j))
    return all(any(set(b) <= s for s in found.values()) for b in blocks)

def attack(name, resfun, n, fixed, cplx, restarts, blocks, kind, scale=1.5):
    nfree = n - len(fixed); dim = (4 if cplx else 2) * nfree
    t0 = time.time(); best = None; good = 0
    for r in range(restarts):
        x0 = rng.normal(0, scale, dim)
        def padded(X):
            v = resfun(X, cplx); return np.concatenate([v, np.zeros(max(0, dim - len(v)))])
        sol = least_squares(padded, x0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        P = unpack(sol.x, n, fixed, cplx); res = np.sqrt(2 * sol.cost)
        sep = min(abs(np.linalg.norm(P[i] - P[j])) for i, j in itertools.combinations(range(n), 2))
        bounded = np.max(np.abs(P)) <= BOUND
        if sep >= SEP and bounded and (best is None or res < best[0]): best = (res, sep, P)
        if sep >= SEP and bounded and res < 1e-10:
            if cplx or structure_ok(P.real, blocks, kind): good += 1
    if best is None:
        print(f"{name} [{'C' if cplx else 'R'}]: restarts={restarts} {time.time()-t0:.0f}s  no run ended with points >= {SEP} apart and |coords| <= {BOUND}"); return
    res, sep, P = best
    print(f"{name} [{'C' if cplx else 'R'}]: restarts={restarts} {time.time()-t0:.0f}s  best residual={res:.2e} (sep {sep:.3f}); "
          f"runs with residual<1e-10, distinct, bounded{'' if cplx else ', structure verified'}: {good}")
    if good and not cplx: print("   !!! REAL REALISATION:", np.round(P.real, 6).tolist())

SQS8 = [(1,2,5,6),(3,4,5,6),(1,3,5,7),(2,4,5,7),(2,3,6,7),(1,4,6,7),(2,3,5,8),(1,4,5,8),(1,3,6,8),(2,4,6,8),(1,2,7,8),(3,4,7,8),(1,2,3,4),(5,6,7,8)]
cube = [frozenset(i - 1 for i in b) for b in SQS8[:12]]
fx3 = {0: (0.0, 0.0), 1: (1.0, 0.0), 2: (0.0, 1.0)}
cand9 = [frozenset(b) for b in [(0,1,2,3,4),(0,5,6,7,8)]] + [frozenset(b) for b in SQS8[:12]]
MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
fano = [frozenset(l) for l in [(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)]]
fxF = {0: (0.0, 0.0), 1: (1.0, 0.0), 3: (0.0, 1.0)}
pappus = [frozenset(l) for l in [(0,1,2),(3,4,5),(6,7,8),(0,4,8),(0,5,7),(1,3,8),(1,5,6),(2,3,7),(2,4,6)]]
R = int(sys.argv[1]) if len(sys.argv) > 1 else 150
attack("control cube n=8 (realisable)", lambda X, c: mob_res(X, cube, fx3, 8, c), 8, fx3, False, 40, cube, "mob")
attack("control Pappus (realisable)", lambda X, c: line_res(X, pappus, fxF, 9, c), 9, fxF, False, 40, pappus, "line")
attack("control Mobius-Kantor over C (realisable over C)", lambda X, c: line_res(X, MK, fx3, 8, c), 8, fx3, True, 40, MK, "line")
attack("control 9-point candidate over C", lambda X, c: mob_res(X, cand9, fx3, 9, c), 9, fx3, True, 40, cand9, "mob")
attack("TARGET Mobius-Kantor (8_3) over R", lambda X, c: line_res(X, MK, fx3, 8, c), 8, fx3, False, R, MK, "line")
attack("TARGET Fano over R", lambda X, c: line_res(X, fano, fxF, 7, c), 7, fxF, False, R, fano, "line")
attack("TARGET 9-point candidate (2x5-block+12 four-blocks) over R", lambda X, c: mob_res(X, cand9, fx3, 9, c), 9, fx3, False, R, cand9, "mob")
sub = [b for b in cand9 if 1 in b]
attack("TARGET only the 7 blocks through point 1 (Fano-derived) over R", lambda X, c: mob_res(X, sub, fx3, 9, c), 9, fx3, False, R, sub, "mob")
sub2 = [b for b in cand9 if 1 not in b]
attack("info: candidate minus the blocks through point 1, over R", lambda X, c: mob_res(X, sub2, fx3, 9, c), 9, fx3, False, 60, sub2, "mob")
