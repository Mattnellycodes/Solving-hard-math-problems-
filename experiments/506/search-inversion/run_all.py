"""Batch driver: subset + inversion-centre optimisation over many rich Moebius universes.

For every universe U and n in NS the CP-SAT model of universe_opt.optimise is solved for every
centre candidate (one representative per point class of U, plus the highest-degree extra sphere
points), and every optimum S is re-evaluated with the full float evaluator of sphere.py (best centre
over ALL pairwise block intersections of S itself), which can only improve the count.

Usage:  python3 run_all.py GROUP [n_min n_max] [time_limit_point time_limit_extra] [max_extras]
        GROUP in {lattice, poly, tri, grid, polygon, pencil, klein, all}
Results are appended to runs/results.jsonl (one line per (universe, n)); runs already present are
skipped (resumable).  Records with count <= f(n) are also appended to runs/records.jsonl.
"""
import sys, os, json, time, math, traceback
import numpy as np
from collections import defaultdict
import sphere as SP
import universes as UV
import lattice_sphere as LS
from universe_opt import Universe, optimise

HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, 'runs')
os.makedirs(RUNS, exist_ok=True)
JOINT = bool(os.environ.get('JOINT'))       # joint subset+centre model (universe_opt2), one solve per n
RESULTS = os.path.join(RUNS, 'results2.jsonl' if JOINT else 'results.jsonl')
RECORDS = os.path.join(RUNS, 'records2.jsonl' if JOINT else 'records.jsonl')
if JOINT:
    from universe_opt2 import optimise_joint


def f(n):
    return SP.formula(n)


# ---------------------------------------------------------------- universe catalogue
def catalogue(group):
    """Yields (name, builder) where builder() -> (name, V (Nx3 unit vectors), labels, exact)."""
    items = []
    if group in ('lattice', 'all'):
        for N in (2, 3, 5, 6, 9, 11, 14, 17, 18, 21, 27, 30, 33, 35, 26, 29, 38):
            items.append((f"S2lat-{N}", lambda N=N: LS.lattice_sphere_universe(N)))
    if group in ('poly', 'all'):
        names = ['cube', 'octahedron', 'cuboctahedron', 'cube+octahedron', 'cube+octa+cubocta',
                 'icosahedron', 'dodecahedron', 'icosidodecahedron', 'icosa+dodeca',
                 'truncated-tetrahedron', 'truncated-cube', 'truncated-octahedron', 'rhombicuboctahedron']
        for m in (3, 4, 5, 6, 8):
            for h in (0.3, 0.5, 1.0, 0.7071):
                names.append(f"prism-{m}-{h}")
                names.append(f"antiprism-{m}-{h}")
        for m in (4, 5, 6, 8, 10):
            names.append(f"bipyramid-{m}")
        for nm in names:
            def b(nm=nm):
                name, P, lab = UV.polyhedron_universe(nm)
                return (name, np.array(P, float), lab, None)
            items.append((nm, b))
    if group in ('tri', 'all'):
        for key in UV.HERONIAN:
            def b(key=key):
                A, B, C = UV.HERONIAN[key]
                from fractions import Fraction as Fr
                A = tuple(Fr(c) for c in A); B = tuple(Fr(c) for c in B); C = tuple(Fr(c) for c in C)
                name, pts, lab, ex = UV.triangle_universe(A, B, C, incentre=True, name=f"tri-{key}")
                return (name, SP.plane_to_sphere(pts), lab, ex)
            items.append((f"tri-{key}", b))
        for i in range(8):
            def b(i=i):
                name, pts, lab, ex = UV.special_triangle_universes()[i]
                return (name, SP.plane_to_sphere(pts), lab, ex)
            items.append((f"tri-special-{i}", b))
    if group in ('grid', 'all'):
        specs = [(4, 4, False), (5, 5, False), (4, 4, True), (3, 3, True), (6, 6, False), (4, 6, False), (3, 5, True), (5, 5, True)]
        for a, bb, half in specs:
            def b(a=a, bb=bb, half=half):
                name, pts, lab, ex = UV.grid_universe(a, bb, half=half, inf=True)
                return (name, SP.plane_to_sphere(pts), lab, ex)
            items.append((f"grid{a}x{bb}{'h' if half else ''}", b))
    if group in ('polygon', 'all'):
        s2, s3, phi = math.sqrt(2), math.sqrt(3), (1 + math.sqrt(5)) / 2
        sec = lambda m: 1 / math.cos(math.pi / m)
        radii = {3: [1, 2, 3, 0.5, s3, 2 * s3, 1.5], 4: [1, 2, 3, s2, 2 * s2, 1 + s2, 4 + math.sqrt(15)],
                 5: [1, 2, phi, phi ** 2, sec(5), sec(5) ** 2, sec(10), 1 / phi], 6: [1, 2, 3, s3, sec(6), 2 * sec(6), 1.5],
                 8: [1, 2, s2, 1 + s2, sec(8), sec(8) ** 2], 10: [1, 2, phi, sec(10), sec(5)], 12: [1, 2, s3, sec(12), sec(6)],
                 7: [1, 2, sec(7), sec(7) ** 2, sec(14)], 9: [1, 2, sec(9), sec(9) ** 2, sec(18)]}
        for m, rs in radii.items():
            def b(m=m, rs=rs):
                rr = []; rots = []
                for r in rs:
                    rr += [r, r]; rots += [0, 1]
                name, pts, lab, ex = UV.polygon_universe(m, rr, rots, centre=True, inf=True, name=f"poly{m}-multi")
                return (name, SP.plane_to_sphere(pts), lab, ex)
            items.append((f"poly{m}-multi", b))
    if group in ('pencil', 'all'):
        specs = [(3.0, 1.0, [0.0, 0.5, -0.5, 1.0], [0.0, math.pi / 4, math.pi / 2, 3 * math.pi / 4]),
                 (2.5, 1.2, [0.0, 0.7, -0.7], [0.0, math.pi / 6, math.pi / 3, math.pi / 2, 2 * math.pi / 3, 5 * math.pi / 6]),
                 (4.0, 2.0, [0.0, 1.0, -1.0, 2.0, -2.0], [math.pi / 2, math.pi / 4, 3 * math.pi / 4])]
        for i, (d, r2, R, ang) in enumerate(specs):
            def b(d=d, r2=r2, R=R, ang=ang, i=i):
                name, pts, lab, ex = UV.two_circle_pencil_universe(d, r2, R, ang, name=f"pencil-{i}")
                return (name, SP.plane_to_sphere(pts), lab, ex)
            items.append((f"pencil-{i}", b))
    if group in ('conic', 'all'):
        specs = [('ellipse24-2-1', lambda: UV.ellipse_universe(24, 2, 1, ('C', 'I'))),
                 ('ellipse24-5-3', lambda: UV.ellipse_universe(24, 5, 3, ('C', 'F', 'I'))),
                 ('ellipse16-2-1', lambda: UV.ellipse_universe(16, 2, 1, ('C', 'I'))),
                 ('parabola13', lambda: UV.parabola_universe(range(-6, 7), ('F', 'I'))),
                 ('parabola13D', lambda: UV.parabola_universe(range(-6, 7), ('F', 'D', 'I'))),
                 ('hyperbola22', lambda: UV.hyperbola_universe(None, ('C', 'I')))]
        for nm, bb in specs:
            def b(bb=bb):
                name, pts, lab, ex = bb()
                return (name, SP.plane_to_sphere(pts), lab, ex)
            items.append((nm, b))
    if group in ('klein', 'all'):
        for m2 in (6, 8, 10, 12):
            for a in (0.3, 0.6):
                def b(m2=m2, a=a):
                    name, pts, lab, ex = UV.klein_polygon_universe(m2, a)
                    return (name, SP.plane_to_sphere(pts), lab, ex)
                items.append((f"klein{m2}-a{a}", b))
    return items


# ---------------------------------------------------------------- extras for big universes
def compute_extras(U, max_pairs=1_500_000, top=40, min_deg=3):
    """Sphere points on >= min_deg blocks that are not universe points, from intersections of pairs
    of rich blocks (size threshold chosen adaptively so that the number of pairs stays bounded)."""
    sizes = np.array([len(b) for b in U.blocks])
    thr = 4
    while True:
        idx = np.flatnonzero(sizes >= thr)
        if len(idx) * (len(idx) - 1) // 2 <= max_pairs or thr >= sizes.max():
            break
        thr += 1
    if len(idx) < 2:
        return np.zeros((0, 3)), np.zeros(0, int), []
    C = SP.candidate_points(U.Nn[idx], U.D[idx], max_pairs=max_pairs)
    if len(C) == 0:
        return np.zeros((0, 3)), np.zeros(0, int), []
    Uc, mult = SP.merge_points(C, min_mult=2)
    if len(Uc) == 0:
        return np.zeros((0, 3)), np.zeros(0, int), []
    dist = np.min(np.linalg.norm(Uc[:, None, :] - U.V[None, :, :], axis=2), axis=1) if len(Uc) < 200000 else None
    if dist is not None:
        Uc = Uc[dist > 1e-6]
    # degrees against all blocks (chunked); cap the number of candidates first by degree on rich blocks
    if len(Uc) > 20000:
        dr = SP.degrees(Uc, U.Nn[idx], U.D[idx])
        keep = np.argsort(-dr)[:20000]
        Uc = Uc[keep]
    deg = SP.degrees(Uc, U.Nn, U.D)
    if dist is None:
        dist = np.min(np.linalg.norm(Uc[:, None, :] - U.V[None, :, :], axis=2), axis=1)
        m = dist > 1e-6; Uc = Uc[m]; deg = deg[m]
    order = np.argsort(-deg)
    order = [int(t) for t in order if deg[t] >= min_deg][:top]
    E = Uc[order]; Ed = deg[order]
    Eb = [SP.blocks_through(x, U.Nn, U.D) for x in E]
    return E, Ed, Eb


def build(name, builder, max_extras):
    nm, V, lab, ex = builder()
    U = Universe(nm, V, lab, ex, extra_centres=False)
    U.extras, U.extra_deg, U.extra_blocks = compute_extras(U, top=max_extras)
    return U


# ---------------------------------------------------------------- full re-evaluation of a subset
def reevaluate(U, S):
    """Best count of the n-subset S over all inversion centres (pairwise block intersections of S's own
    blocks, plus universe points), float evaluation.  Returns (count, centre_xyz, deg)."""
    V = U.V[S]
    ev = SP.evaluate(V)
    if ev is None:
        return None
    best = ev.best_count; deg = ev.best_deg
    o = ev.best_centres(1)[0][1] if len(ev.cdeg) else None
    # universe points as centres
    Sset = set(S)
    for p in range(U.N):
        if p in Sset:
            continue
        d = len(SP.blocks_through(U.V[p], ev.Nn, ev.D))
        if ev.nb - d < best:
            best = ev.nb - d; deg = d; o = U.V[p]
    return best, (None if o is None else [float(c) for c in o]), deg, ev.nb


def done_keys():
    keys = set()
    if os.path.exists(RESULTS):
        for line in open(RESULTS):
            try:
                r = json.loads(line); keys.add((r['universe'], r['n']))
            except Exception:
                pass
    return keys


def run(group, ns, tl_point, tl_extra, max_extras):
    done = done_keys()
    only = os.environ.get('UNIV')          # optional: run a single universe by name
    for name, builder in catalogue(group):
        if only and name != only:
            continue
        if all((name, n) in done for n in ns):
            print(f"[{name}] all done, skipping", flush=True); continue
        t0 = time.time()
        try:
            U = build(name, builder, max_extras)
        except Exception as e:
            print(f"[{name}] BUILD FAILED: {e}", flush=True); traceback.print_exc(); continue
        info = U.info()
        classes = U.point_classes()
        reps = sorted([cl[0] for cl in classes], key=lambda i: -U.pdeg[i])
        print(f"[{name}] {info} classes={len(classes)} extras={len(U.extras)} maxextra={info['max_extra_deg']} build={time.time()-t0:.1f}s", flush=True)
        centres = [('point', i) for i in reps] + [('extra', j) for j in range(len(U.extras))]
        for n in ns:
            if (name, n) in done or n + 1 > U.N:
                continue
            t1 = time.time(); best = None
            if JOINT:
                try:
                    st, cnt, S, c, bound = optimise_joint(U, n, time_limit=tl_point, workers=2)
                except Exception as e:
                    print(f"   optimise_joint failed: {e}", flush=True); cnt = None
                if cnt is not None:
                    best = dict(count=cnt, status=st, bound=bound, centre=list(c), S=S)
            else:
                for c in centres:
                    tl = tl_point if c[0] == 'point' else tl_extra
                    try:
                        st, cnt, S, bound = optimise(U, n, c, time_limit=tl, workers=2)
                    except Exception as e:
                        print(f"   optimise failed {c}: {e}", flush=True); continue
                    if cnt is None:
                        continue
                    rec = dict(count=cnt, status=st, bound=bound, centre=list(c), S=S)
                    if best is None or cnt < best['count']:
                        best = rec
            if best is None:
                continue
            # full re-evaluation of the best subset (and of the runner-up centres' subsets is skipped)
            re = reevaluate(U, best['S'])
            out = dict(universe=name, n=n, N=U.N, formula=f(n), cp_count=best['count'], cp_status=best['status'],
                       centre=best['centre'], centre_label=(U.labels[best['centre'][1]] if best['centre'][0] == 'point'
                                                             else f"extra{best['centre'][1]}(deg {int(U.extra_deg[best['centre'][1]])})"),
                       S=best['S'], S_labels=[U.labels[p] for p in best['S']], seconds=round(time.time() - t1, 1))
            if re is not None:
                out.update(count=re[0], best_centre_xyz=re[1], best_deg=re[2], nblocks_S=re[3])
            else:
                out.update(count=best['count'])
            flag = " <<< BELOW FORMULA" if out['count'] < f(n) else (" (= formula)" if out['count'] == f(n) else "")
            print(f"   n={n}: count={out['count']} (cp {best['count']}, f={f(n)}) centre={out['centre_label']} S={out['S_labels']} t={out['seconds']}s{flag}", flush=True)
            with open(RESULTS, 'a') as fh:
                fh.write(json.dumps(out) + "\n")
            if out['count'] <= f(n):
                out2 = dict(out); out2['S_xyz'] = [[float(c) for c in U.V[p]] for p in best['S']]
                if U.exact is not None:
                    ex = U.exact
                    if 'points3d' in ex:
                        out2['exact'] = {'points3d': [ex['points3d'][p] for p in best['S']], 'N': ex['N']}
                    else:
                        out2['exact'] = {'planar': [str(ex['planar'][p]) for p in best['S']], 'radicands': list(ex['radicands'])}
                with open(RECORDS, 'a') as fh:
                    fh.write(json.dumps(out2) + "\n")


if __name__ == "__main__":
    group = sys.argv[1] if len(sys.argv) > 1 else 'all'
    n0 = int(sys.argv[2]) if len(sys.argv) > 2 else 9
    n1 = int(sys.argv[3]) if len(sys.argv) > 3 else 16
    tlp = float(sys.argv[4]) if len(sys.argv) > 4 else 20.0
    tle = float(sys.argv[5]) if len(sys.argv) > 5 else 8.0
    mx = int(sys.argv[6]) if len(sys.argv) > 6 else 30
    run(group, list(range(n0, n1 + 1)), tlp, tle, mx)
