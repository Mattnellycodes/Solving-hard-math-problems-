"""Numerical realisation attack, second (clean) version.  Levenberg-Marquardt on the incidence
equations only (no barrier); distinctness is checked afterwards.  Random restarts.
Real and COMPLEX coordinates: the complex runs are positive controls showing that the solver finds
the (known) complex realisations of Mobius-Kantor and of the Fano-type structures, so a failure of the
real runs is not a solver artefact.
Mobius structures: block = concyclic-or-collinear quadruples: det[x^2+y^2, x, y, 1] = 0 for every
4-subset (complex case: same polynomial, x,y complex).  Three points fixed at (0,0),(1,0),(0,1).
Line structures: det[x,y,1]=0 per triple; four frame points fixed (no three on a structure line)."""
import numpy as np, itertools, sys, time
from scipy.optimize import least_squares
rng = np.random.default_rng(2024)
SEP = 0.05   # minimum pairwise point distance for a solution to count as a realisation with distinct points

def unpack(X, n, fixed, cplx):
    P = np.zeros((n, 2), dtype=complex if cplx else float)
    for i, v in fixed.items(): P[i] = v
    free = [i for i in range(n) if i not in fixed]
    if cplx:
        X = X[: 2 * len(free)] + 1j * X[2 * len(free):]
    P[free] = X.reshape(-1, 2)
    return P

def mob_res(X, blocks, fixed, n, cplx):
    P = unpack(X, n, fixed, cplx)
    out = []
    for b in blocks:
        for q in itertools.combinations(sorted(b), 4):
            M = np.array([[P[i, 0] ** 2 + P[i, 1] ** 2, P[i, 0], P[i, 1], 1.0] for i in q])
            out.append(np.linalg.det(M))
    out = np.array(out)
    return np.concatenate([out.real, out.imag]) if cplx else out

def line_res(X, lines, fixed, m, cplx):
    P = unpack(X, m, fixed, cplx)
    out = []
    for l in lines:
        for t in itertools.combinations(sorted(l), 3):
            M = np.array([[P[i, 0], P[i, 1], 1.0] for i in t])
            out.append(np.linalg.det(M))
    out = np.array(out)
    return np.concatenate([out.real, out.imag]) if cplx else out

def min_sep(P):
    return min(abs(np.linalg.norm(P[i] - P[j])) for i, j in itertools.combinations(range(len(P)), 2))

def attack(name, resfun, n, fixed, cplx, restarts, scale=2.0):
    nfree = n - len(fixed); dim = (4 if cplx else 2) * nfree
    best = None; t0 = time.time(); good = 0
    for r in range(restarts):
        x0 = rng.normal(0, scale, dim)
        def padded(X):
            r = resfun(X, cplx)
            return np.concatenate([r, np.zeros(max(0, dim - len(r)))])
        sol = least_squares(padded, x0, method="lm", xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=3000)
        P = unpack(sol.x, n, fixed, cplx)
        sep = min_sep(P); res = np.sqrt(2 * sol.cost)
        # scale-normalised residual: divide by (1+max|P|)^4 to be fair to large coordinates
        nres = res / (1 + np.max(np.abs(P))) ** 4
        if sep > SEP and (best is None or nres < best[0]): best = (nres, res, sep, P)
        if sep > SEP and nres < 1e-10: good += 1
    if best is None:
        print(f"{name} [{'complex' if cplx else 'real'}]: restarts={restarts} time={time.time()-t0:.1f}s  NO solution with all points >= {SEP} apart at all (every restart collapsed points)")
        return None
    nres, res, sep, P = best
    print(f"{name} [{'complex' if cplx else 'real'}]: restarts={restarts} time={time.time()-t0:.1f}s  "
          f"best normalised residual={nres:.2e} (raw {res:.2e}), min point separation={sep:.3f}, "
          f"#restarts with residual<1e-10 and all points >= {SEP} apart: {good}")
    if good and not cplx:
        print("   !!! REAL REALISATION FOUND:", np.round(P, 6).tolist())
    return best

SQS8 = [(1,2,5,6),(3,4,5,6),(1,3,5,7),(2,4,5,7),(2,3,6,7),(1,4,6,7),(2,3,5,8),(1,4,5,8),
        (1,3,6,8),(2,4,6,8),(1,2,7,8),(3,4,7,8),(1,2,3,4),(5,6,7,8)]
cube = [frozenset(i - 1 for i in b) for b in SQS8[:12]]
fx3 = {0: (0.0, 0.0), 1: (1.0, 0.0), 2: (0.0, 1.0)}
cand9 = [frozenset(b) for b in [(0,1,2,3,4),(0,5,6,7,8)]] + [frozenset(i for i in b) for b in SQS8[:12]]
# sanity: cand9 is the enumeration's unique survivor structure (two 5-blocks sharing point 0 + 12 four-blocks
# on {1..8}); check pairwise intersections <= 2 and degrees
assert all(len(a & b) <= 2 for a, b in itertools.combinations(cand9, 2))
print("cand9 degrees:", [sum(1 for b in cand9 if p in b) for p in range(9)])
MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
fxMK = {0: (0.0, 0.0), 1: (1.0, 0.0), 2: (0.0, 1.0)}   # {0,1,2} is not an MK line
fano = [frozenset(l) for l in [(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)]]
fxF = {0: (0.0, 0.0), 1: (1.0, 0.0), 3: (0.0, 1.0)}
assert not any(len(set(fxF) & l) >= 3 for l in fano)
pappus = [frozenset(l) for l in [(0,1,2),(3,4,5),(6,7,8),(0,4,8),(0,5,7),(1,3,8),(1,5,6),(2,3,7),(2,4,6)]]
fxP = {0: (0.0, 0.0), 1: (1.0, 0.0), 3: (0.0, 1.0)}

which = sys.argv[1] if len(sys.argv) > 1 else "all"
R = int(sys.argv[2]) if len(sys.argv) > 2 else 150
if which in ("controls", "all"):
    attack("control: cube structure n=8 (realisable)", lambda X, c: mob_res(X, cube, fx3, 8, c), 8, fx3, False, 40)
    attack("control: Pappus 9_3 (realisable)", lambda X, c: line_res(X, pappus, fxP, 9, c), 9, fxP, False, 40)
    attack("control: Mobius-Kantor over C (realisable over C)", lambda X, c: line_res(X, MK, fxMK, 8, c), 8, fxMK, True, 40)
    attack("control: Fano over C? (no: Fano needs char 2) ", lambda X, c: line_res(X, fano, fxF, 7, c), 7, fxF, True, 40)
if which in ("targets", "all"):
    attack("target: Mobius-Kantor (8_3) lines", lambda X, c: line_res(X, MK, fxMK, 8, c), 8, fxMK, False, R)
    attack("target: Fano plane lines", lambda X, c: line_res(X, fano, fxF, 7, c), 7, fxF, False, R)
    attack("target: 9-point candidate structure (2x5-block + 12 four-blocks), real", lambda X, c: mob_res(X, cand9, fx3, 9, c), 9, fx3, False, R)
    attack("control: 9-point candidate structure over C", lambda X, c: mob_res(X, cand9, fx3, 9, c), 9, fx3, True, 40)
    sub = [b for b in cand9 if 1 in b]
    attack("target: the 7 blocks through point 1 only (Fano-derived), real", lambda X, c: mob_res(X, sub, fx3, 9, c), 9, fx3, False, R)
    sub2 = [b for b in cand9 if 1 not in b]
    attack("info: candidate minus blocks through point 1, real", lambda X, c: mob_res(X, sub2, fx3, 9, c), 9, fx3, False, 60)
