"""1-parameter family scanner for Erdos #506 (sphere model, floats, vectorised).

A family t -> P(t) (planar points, 'inf' allowed) has a fixed block structure except at finitely
many "events":
  (E1) four points become co-circular (a new block / a bigger block): sign change of the 4x4
       coplanarity determinant of their sphere images;
  (E2) three blocks become concurrent at a sphere point that is not a configuration point: sign
       change of the signed distance from an intersection point of two blocks to a third block.
The count of the best planar picture of P(t) is |B(t)| - max_O deg(O), so it can only improve at
events of type (E1) (more/bigger blocks) or (E2) (a centre of higher degree).  We scan t on a grid,
detect sign changes, refine them by bisection and evaluate the exact-tolerance count at the event.
"""
import math, itertools, sys, json, time
import numpy as np
import sphere as SP


def state(V, tol=SP.TOL_ANGLE):
    """Block keys, the intersection points of all block pairs (with a canonical sign), and the
    residual matrix R[(pair, sign), block]."""
    blocks = SP.find_blocks(V, tol)
    keys = tuple(tuple(sorted(b)) for b in blocks)
    Nn, D = SP.planes(V, blocks)
    B = len(blocks)
    if B < 3:
        return keys, None, None
    I, J = np.triu_indices(B, 1)
    n1 = Nn[I]; n2 = Nn[J]; d1 = D[I]; d2 = D[J]
    u = np.cross(n1, n2); un = np.linalg.norm(u, axis=1)
    ok = un > 1e-9
    I = I[ok]; J = J[ok]; n1 = n1[ok]; n2 = n2[ok]; d1 = d1[ok]; d2 = d2[ok]; u = u[ok] / un[ok][:, None]
    # canonical orientation of u
    sgn = np.ones(len(u))
    first = np.argmax(np.abs(u) > 1e-9, axis=1)
    sgn = np.sign(u[np.arange(len(u)), first]); sgn[sgn == 0] = 1
    u = u * sgn[:, None]
    g = np.sum(n1 * n2, axis=1); det = 1 - g * g
    a = (d1 - d2 * g) / det; b = (d2 - d1 * g) / det
    x0 = a[:, None] * n1 + b[:, None] * n2
    h2 = 1 - np.sum(x0 * x0, axis=1)
    ok2 = h2 > -1e-9
    I = I[ok2]; J = J[ok2]; x0 = x0[ok2]; u = u[ok2]; h = np.sqrt(np.maximum(h2[ok2], 0))
    X = np.vstack([x0 + h[:, None] * u, x0 - h[:, None] * u])
    R = X @ Nn.T - D[None, :]
    # mask: own blocks, and configuration points
    P = len(I)
    rows = np.arange(2 * P); pi = np.concatenate([I, I]); pj = np.concatenate([J, J])
    R[rows, pi] = np.nan; R[rows, pj] = np.nan
    dist = np.min(np.linalg.norm(X[:, None, :] - V[None, :, :], axis=2), axis=1)
    R[dist < 1e-6, :] = np.nan
    return keys, (I, J), R


def det4(V):
    N = len(V)
    M = np.hstack([V, np.ones((N, 1))])
    Q = np.array(list(itertools.combinations(range(N), 4)))
    return Q, np.linalg.det(M[Q])


def scan(family, t_grid, name="family", verbose=True, max_events=400, eval_extra=None):
    """family(t) -> planar point list.  Returns (events, best)."""
    events = []; seen_t = set()
    prev = None
    n = len(family(t_grid[0]))
    fval = SP.formula(n)
    best = (10 ** 9, None)
    Q4 = None

    def evaluate_event(ts, kind, key):
        nonlocal best
        tk = round(ts, 9)
        if tk in seen_t:
            return
        seen_t.add(tk)
        ev = SP.evaluate_planar(family(ts))
        if ev is None:
            return
        s = ev.summary()
        e = dict(t=ts, kind=kind, key=str(key), count=s['best_count'], nblocks=s['nblocks'], sizes=s['sizes'], best_deg=s['best_deg'])
        events.append(e)
        if s['best_count'] < best[0]:
            best = (s['best_count'], e)
            if verbose:
                flag = " <<< BELOW FORMULA" if s['best_count'] < fval else (" (= formula)" if s['best_count'] == fval else "")
                print(f"[{name}] t={ts:.12f} {kind} count={s['best_count']} (f={fval}) blocks={s['nblocks']} sizes={s['sizes']} deg={s['best_deg']}{flag}", flush=True)

    for t in t_grid:
        V = SP.plane_to_sphere(family(t))
        Q4, dets = det4(V)
        keys, pairs, R = state(V)
        cur = (t, dets, keys, pairs, R)
        if prev is not None:
            t0, dets0, keys0, pairs0, R0 = prev
            # E1: coplanarity sign changes
            ch = np.flatnonzero((dets0 * dets < 0) & (np.abs(dets0) > 1e-13) & (np.abs(dets) > 1e-13))
            for c in ch[:max_events]:
                q = tuple(Q4[c])
                lo, hi = t0, t
                def val(tt):
                    Vt = SP.plane_to_sphere(family(tt))
                    M = np.hstack([Vt, np.ones((len(Vt), 1))])
                    return float(np.linalg.det(M[list(q)]))
                vlo = val(lo)
                for _ in range(50):
                    mid = 0.5 * (lo + hi); vm = val(mid)
                    if vlo * vm <= 0:
                        hi = mid
                    else:
                        lo, vlo = mid, vm
                    if hi - lo < 1e-13 * max(1.0, abs(hi)):
                        break
                evaluate_event(0.5 * (lo + hi), 'E1', q)
            # E2: concurrency sign changes (only if the block structure is unchanged)
            if keys == keys0 and R is not None and R0 is not None and R.shape == R0.shape:
                M = (R0 * R < 0) & (np.abs(R0) > 1e-12) & (np.abs(R) > 1e-12)
                rows, cols = np.nonzero(M)
                found = []          # (t*, residual matrix at t*) of events already located in this interval
                for r_, c_ in list(zip(rows, cols))[:max_events]:
                    # skip if this residual already vanishes at an event found in this interval
                    if any(abs(Rf[r_, c_]) < 1e-7 for tf, Rf in found if Rf is not None):
                        continue
                    lo, hi = t0, t
                    def val(tt):
                        Vt = SP.plane_to_sphere(family(tt))
                        k2, p2, R2 = state(Vt)
                        if k2 != keys or R2 is None or R2.shape != R.shape:
                            return float('nan')
                        return float(R2[r_, c_])
                    vlo = val(lo); bad = False
                    for _ in range(50):
                        mid = 0.5 * (lo + hi); vm = val(mid)
                        if vm != vm:
                            bad = True; break
                        if vlo * vm <= 0:
                            hi = mid
                        else:
                            lo, vlo = mid, vm
                        if hi - lo < 1e-13 * max(1.0, abs(hi)):
                            break
                    if bad:
                        continue
                    ts = 0.5 * (lo + hi)
                    k2, p2, R2 = state(SP.plane_to_sphere(family(ts)))
                    found.append((ts, R2 if (k2 == keys and R2 is not None and R2.shape == R.shape) else None))
                    evaluate_event(ts, 'E2', (int(r_), int(c_)))
        prev = cur
    return events, best


# ---------------------------------------------------------------- families
def two_polygons(m, rot, extra=()):
    """Two concentric regular m-gons, radii 1 and t, second rotated by rot*pi/m; extra in {'C','I'}."""
    def fam(t):
        pts = []
        for k in range(m):
            th = 2 * math.pi * k / m
            pts.append((math.cos(th), math.sin(th)))
        for k in range(m):
            th = 2 * math.pi * k / m + rot * math.pi / m
            pts.append((t * math.cos(th), t * math.sin(th)))
        if 'C' in extra:
            pts.append((0.0, 0.0))
        if 'I' in extra:
            pts.append('inf')
        return pts
    return fam


def orthocentric_family(a, extra=('feet', 'O', 'N', 'mid', 'inf')):
    """Triangle A=(0,a), B=(-1,0), C=(t,0): orthocentre H, feet, circumcentre O, nine-point centre N,
    midpoints, Euler points, reflections of H, centroid, plus infinity (1-parameter slice in t)."""
    def fam(t):
        A = np.array([0.0, a]); B = np.array([-1.0, 0.0]); C = np.array([t, 0.0])
        def circum(P, Q, R):
            ax, ay = P; bx, by = Q; cx, cy = R
            d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
            a2 = ax * ax + ay * ay; b2 = bx * bx + by * by; c2 = cx * cx + cy * cy
            return np.array([(a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d, (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d])
        O = circum(A, B, C); H = A + B + C - 2 * O
        pts = [A, B, C, H]
        def foot(P, Q, R):
            d = R - Q; return Q + d * ((P - Q) @ d) / (d @ d)
        if 'feet' in extra:
            pts += [foot(A, B, C), foot(B, C, A), foot(C, A, B)]
        if 'O' in extra:
            pts.append(O)
        if 'N' in extra:
            pts.append((O + H) / 2)
        if 'mid' in extra:
            pts += [(B + C) / 2, (C + A) / 2, (A + B) / 2]
        if 'euler' in extra:
            pts += [(A + H) / 2, (B + H) / 2, (C + H) / 2]
        if 'reflH' in extra:
            pts += [2 * foot(A, B, C) - H, 2 * foot(B, C, A) - H, 2 * foot(C, A, B) - H]
        if 'G' in extra:
            pts.append((A + B + C) / 3)
        out = [(float(p[0]), float(p[1])) for p in pts]
        if 'inf' in extra:
            out.append('inf')
        ded = []
        for p in out:
            if p == 'inf' or all(q == 'inf' or abs(p[0] - q[0]) > 1e-9 or abs(p[1] - q[1]) > 1e-9 for q in ded):
                ded.append(p)
        return ded
    return fam


def box_family(b, extra=()):
    """Rectangular box 1 x b x t inscribed in a sphere (cube structure) + optional 'faces' (the 6
    axis points of the sphere); returned as a generic stereographic picture."""
    def fam(t):
        P = [(sx, sy * b, sz * t) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
        R = math.sqrt(1 + b * b + t * t)
        if 'faces' in extra:
            P += [(R, 0, 0), (-R, 0, 0), (0, R, 0), (0, -R, 0), (0, 0, R), (0, 0, -R)]
        V = np.array(P, float) / R
        o = np.array([0.31, 0.47, 0.83]); o /= np.linalg.norm(o)
        pl = SP.sphere_to_plane(V, o)
        return [(float(x), float(y)) for x, y in pl]
    return fam


def run_family(fam, grid, name, out):
    t0 = time.time()
    ev, best = scan(fam, grid, name=name)
    n = len(fam(grid[len(grid) // 2]))
    rec = dict(family=name, n=n, formula=SP.formula(n), best=best[0], best_event=best[1], nevents=len(ev), seconds=round(time.time() - t0, 1))
    out.append(rec)
    print("   ==>", {k: v for k, v in rec.items() if k != 'best_event'}, flush=True)
    return rec


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else 'pent'
    out = []
    if which == 'pent':
        for rot in (0, 1):
            for extra in ((), ('C',), ('I',), ('C', 'I')):
                run_family(two_polygons(5, rot, extra), np.linspace(0.02, 0.995, 1500), f"2x5gon rot={rot} extra={extra}", out)
    elif which == 'polys':
        for m in (6, 7, 8):
            for rot in (0, 1):
                for extra in ((), ('C',), ('I',), ('C', 'I')):
                    run_family(two_polygons(m, rot, extra), np.linspace(0.02, 0.995, 1200), f"2x{m}gon rot={rot} extra={extra}", out)
    elif which == 'ortho':
        for a in (1.0, 1.5, 2.0, 3.0, 0.75):
            for extra in (('feet', 'O', 'inf'), ('feet', 'N', 'inf'), ('feet', 'O', 'N', 'inf'), ('feet', 'mid', 'inf'),
                          ('feet', 'euler', 'inf'), ('feet', 'reflH', 'inf'), ('feet', 'O', 'N', 'G', 'inf'), ('feet', 'mid', 'euler', 'inf'),
                          ('feet', 'O', 'reflH', 'inf')):
                run_family(orthocentric_family(a, extra), np.linspace(0.05, 6.0, 1200), f"ortho a={a} extra={extra}", out)
    elif which == 'box':
        for b in (1.0, 1.5, 2.0, math.sqrt(2), math.sqrt(3)):
            for extra in ((), ('faces',)):
                run_family(box_family(b, extra), np.linspace(0.1, 4.0, 1200), f"box b={b} extra={extra}", out)
    with open(f"runs/family_{which}.json", 'w') as fh:
        json.dump(out, fh, indent=1, default=str)
