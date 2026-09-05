"""Codimension-2 events in 2-parameter families (Erdos #506).

Given a family (a, t) -> P(a, t) and a codimension-1 event located at (a0, t0) by family_scan (an E1
coplanarity key = 4 point indices, or an E2 concurrency key = (block_i, block_j, sign, block_k) with
blocks as sorted point tuples), we follow the event curve t*(a) for a on a grid (bisection on the
event's residual, bracket around the previous t*), evaluate the full count at every (a, t*(a)), and
detect secondary events along the curve (sign changes of any other coplanarity determinant or
concurrency residual between consecutive a-values), which are refined by bisection in a.  A secondary
event is a codimension-2 point (a*, t*) where two coincidences hold simultaneously; the count there
can only be lower than on the curve.
"""
import math, itertools, sys, json, time
import numpy as np
import sphere as SP
import family_scan as FS


def residual_E2(V, key):
    """Residual of the E2 key (ki, kj, s, kk) at the sphere set V, or nan if the blocks are absent."""
    keys, pairs, R = FS.state(V)
    if R is None:
        return float('nan'), keys, pairs, R
    ki, kj, s, kk = key
    try:
        i = keys.index(ki); j = keys.index(kj); k = keys.index(kk)
    except ValueError:
        return float('nan'), keys, pairs, R
    I, J = pairs
    if i > j:
        i, j = j, i
    m = np.flatnonzero((I == i) & (J == j))
    if len(m) == 0:
        return float('nan'), keys, pairs, R
    row = int(m[0]) + (0 if s == 1 else len(I))
    return float(R[row, k]), keys, pairs, R


def residual(family2, a, t, kind, key):
    V = SP.plane_to_sphere(family2(a, t))
    if kind == 'E1':
        M = np.hstack([V, np.ones((len(V), 1))])
        return float(np.linalg.det(M[list(key)]))
    return residual_E2(V, key)[0]


def locate_t(family2, a, kind, key, t_lo, t_hi):
    """Bisection for the event's t at parameter a inside [t_lo, t_hi]; returns t* or None."""
    vlo = residual(family2, a, t_lo, kind, key); vhi = residual(family2, a, t_hi, kind, key)
    if vlo != vlo or vhi != vhi or vlo * vhi > 0:
        return None
    for _ in range(60):
        mid = 0.5 * (t_lo + t_hi); vm = residual(family2, a, mid, kind, key)
        if vm != vm:
            return None
        if vlo * vm <= 0:
            t_hi, vhi = mid, vm
        else:
            t_lo, vlo = mid, vm
        if t_hi - t_lo < 1e-13 * max(1.0, abs(t_hi)):
            break
    return 0.5 * (t_lo + t_hi)


def e2_key_from_indices(keys, pairs, row, col):
    I, J = pairs
    P = len(I)
    s = 1 if row < P else -1
    r = row if row < P else row - P
    return (keys[int(I[r])], keys[int(J[r])], s, keys[int(col)])


def follow(family2, kind, key, a0, t0, a_grid, dt=0.05, name="follow", verbose=True, max_secondary=60):
    """Follow the codim-1 event from (a0, t0) over a_grid (which should contain a0 or start near it)."""
    out = dict(name=name, kind=kind, key=str(key), curve=[], secondary=[], best=None)
    n = len(family2(a0, t0)); fval = SP.formula(n)
    best = (10 ** 9, None)
    prev = None
    # order the grid so that we start at the a closest to a0 and expand outwards
    a_grid = list(a_grid)
    a_grid.sort(key=lambda a: abs(a - a0))
    # process in two sweeps (a >= a0 increasing, a < a0 decreasing) with continuity
    for sweep in ([a for a in sorted(a_grid) if a >= a0], [a for a in sorted(a_grid, reverse=True) if a < a0]):
        t_prev = t0; prev = None
        for a in sweep:
            width = dt
            ts = None
            for _ in range(6):
                ts = locate_t(family2, a, kind, key, t_prev - width, t_prev + width)
                if ts is not None:
                    break
                width *= 2
            if ts is None:
                if verbose:
                    print(f"[{name}] lost the event at a={a:.6f}", flush=True)
                break
            t_prev = ts
            V = SP.plane_to_sphere(family2(a, ts))
            Q4, dets = FS.det4(V)
            keys, pairs, R = FS.state(V)
            ev = SP.evaluate(V); s = ev.summary() if ev is not None else None
            cnt = s['best_count'] if s else None
            out['curve'].append(dict(a=a, t=ts, count=cnt, deg=(s['best_deg'] if s else None), nblocks=(s['nblocks'] if s else None)))
            if s and cnt < best[0]:
                best = (cnt, dict(a=a, t=ts, count=cnt, sizes=s['sizes'], deg=s['best_deg']))
                if verbose:
                    flag = " <<< BELOW FORMULA" if cnt < fval else (" (= formula)" if cnt == fval else "")
                    print(f"[{name}] on curve a={a:.9f} t={ts:.9f} count={cnt} (f={fval}) sizes={s['sizes']} deg={s['best_deg']}{flag}", flush=True)
            # secondary events between prev and current
            if prev is not None:
                a_p, t_p, dets_p, keys_p, pairs_p, R_p = prev
                sec = []
                ch = np.flatnonzero((dets_p * dets < 0) & (np.abs(dets_p) > 1e-13) & (np.abs(dets) > 1e-13))
                for c in ch[:max_secondary]:
                    q = tuple(int(v) for v in Q4[c])
                    if kind == 'E1' and q == tuple(key):
                        continue
                    sec.append(('E1', q))
                if keys == keys_p and R is not None and R_p is not None and R.shape == R_p.shape:
                    M = (R_p * R < 0) & (np.abs(R_p) > 1e-12) & (np.abs(R) > 1e-12)
                    rows, cols = np.nonzero(M)
                    for r_, c_ in list(zip(rows, cols))[:max_secondary]:
                        k2 = e2_key_from_indices(keys, pairs, int(r_), int(c_))
                        sec.append(('E2', k2))
                found_here = []
                for kind2, key2 in sec:
                    # skip if already vanishing at a found secondary point
                    if any(abs(residual(family2, af, tf, kind2, key2)) < 1e-8 for af, tf in found_here):
                        continue
                    # bisection in a; at each a locate t on the primary curve
                    lo, hi = (a_p, a) if a_p < a else (a, a_p)
                    tlo, thi = (t_p, ts) if a_p < a else (ts, t_p)
                    def sec_val(aa, tguess):
                        tt = None; w = dt
                        for _ in range(6):
                            tt = locate_t(family2, aa, kind, key, tguess - w, tguess + w)
                            if tt is not None:
                                break
                            w *= 2
                        if tt is None:
                            return float('nan'), None
                        return residual(family2, aa, tt, kind2, key2), tt
                    vlo, _ = sec_val(lo, tlo); vhi, _ = sec_val(hi, thi)
                    if vlo != vlo or vhi != vhi or vlo * vhi > 0:
                        continue
                    ok = True; tmid = None
                    for _ in range(50):
                        mid = 0.5 * (lo + hi); vm, tmid = sec_val(mid, 0.5 * (tlo + thi))
                        if vm != vm:
                            ok = False; break
                        if vlo * vm <= 0:
                            hi, vhi, thi = mid, vm, tmid
                        else:
                            lo, vlo, tlo = mid, vm, tmid
                        if hi - lo < 1e-12 * max(1.0, abs(hi)):
                            break
                    if not ok or tmid is None:
                        continue
                    a2 = 0.5 * (lo + hi); t2 = 0.5 * (tlo + thi)
                    found_here.append((a2, t2))
                    ev2 = SP.evaluate_planar(family2(a2, t2))
                    if ev2 is None:
                        continue
                    s2 = ev2.summary()
                    rec = dict(a=a2, t=t2, kind=kind2, key=str(key2), count=s2['best_count'], sizes=s2['sizes'], deg=s2['best_deg'], nblocks=s2['nblocks'])
                    out['secondary'].append(rec)
                    if s2['best_count'] < best[0]:
                        best = (s2['best_count'], rec)
                        if verbose:
                            flag = " <<< BELOW FORMULA" if s2['best_count'] < fval else (" (= formula)" if s2['best_count'] == fval else "")
                            print(f"[{name}] SECONDARY {kind2} a={a2:.9f} t={t2:.9f} count={s2['best_count']} (f={fval}) sizes={s2['sizes']} deg={s2['best_deg']}{flag}", flush=True)
            prev = (a, ts, dets, keys, pairs, R)
    out['best'] = best[1]
    return out


# ---------------------------------------------------------------- 2-parameter families
def ortho2(extra):
    def fam(a, t):
        return FS.orthocentric_family(a, extra)(t)
    return fam


def polygons2(m, rot, extra=()):
    """two concentric m-gons: radius ratio t, second rotated by a (radians)."""
    def fam(a, t):
        pts = []
        for k in range(m):
            th = 2 * math.pi * k / m
            pts.append((math.cos(th), math.sin(th)))
        for k in range(m):
            th = 2 * math.pi * k / m + a
            pts.append((t * math.cos(th), t * math.sin(th)))
        if 'C' in extra:
            pts.append((0.0, 0.0))
        if 'I' in extra:
            pts.append('inf')
        return pts
    return fam


def box2(extra=()):
    def fam(a, t):
        return FS.box_family(a, extra)(t)
    return fam


if __name__ == "__main__":
    # test: the box family 1 x a x t; the E2 event of the n=8 record at (a=1, t=sqrt3) traces the
    # curve t = sqrt(3) a? (c^2 = 3 b^2) -- follow it and look for secondary events.
    fam = box2(())
    a0 = 1.0
    # find the event at a0 with the 1-D scanner
    evs, best = FS.scan(lambda t: fam(a0, t), np.linspace(1.5, 2.0, 60), name="box slice", verbose=False)
    e = min(evs, key=lambda e: e['count'])
    print("slice event:", e)
    # reconstruct the E2 key from the (row, col) indices at the event
    V = SP.plane_to_sphere(fam(a0, e['t'] + 1e-4))
    keys, pairs, R = FS.state(V)
    row, col = json.loads(e['key'].replace('(', '[').replace(')', ']'))
    key = e2_key_from_indices(keys, pairs, row, col)
    res = follow(fam, 'E2', key, a0, e['t'], np.linspace(0.5, 2.5, 80), dt=0.05, name="box follow")
    print("best:", res['best'], "secondary events:", len(res['secondary']))
