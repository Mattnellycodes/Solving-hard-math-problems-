"""Construction attack v2 (fixes v1's artefacts): closure pools (line/circle intersections of random
rational seeds) or circle pools (rational unit-circle points + centre + chord intersections), simulated
annealing over 9-subsets minimising the float circle count.  Fixes: (1) 9-sets with two points closer
than 2% of their diameter are rejected (near-coincident points fooled the tolerance clustering);
(2) every hit with <= 24 circles is re-verified: exactly (Fractions, own exact counter mycount.py)
when all its points are exact rationals, otherwise by a recount with a 1000x tighter tolerance."""
import numpy as np, itertools, math, sys, time, random
from fractions import Fraction as Fr
from closure_search import block_vectors, TRIPLES, rand_rational
from mycount import analyse

def count_circles(P, tol=1e-6):
    d = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)
    diam = d.max(); np.fill_diagonal(d, np.inf)
    if d.min() < 0.02 * diam: return None
    v, nrm = block_vectors(P)
    if np.any(nrm < 1e-12): return None
    keys = {}
    for i in range(84):
        k = tuple(np.round(v[i] / tol).astype(np.int64))
        keys.setdefault(k, set()).update(TRIPLES[i].tolist())
    if any(len(s) == 9 for s in keys.values()): return None
    circles = sum(1 for k in keys if abs(k[0]) > 0)
    return circles, len(keys) - circles

def exact_line_line(l1, l2):
    (A1, B1, C1), (A2, B2, C2) = l1, l2
    det = A1 * B2 - A2 * B1
    if det == 0: return None
    return ((C1 * B2 - C2 * B1) / det, (A1 * C2 - A2 * C1) / det)

def seed_pool(nseed, nmax, rng_seed):
    """pool entries: (x, y, exact) with exact = (Fr,Fr) or None"""
    random.seed(rng_seed)
    pts = []
    while len(pts) < nseed:
        p = (rand_rational(), rand_rational())
        if p not in pts: pts.append(p)
    pool = [(float(x), float(y), (x, y)) for x, y in pts]
    def add(x, y, ex=None):
        if abs(x) > 40 or abs(y) > 40: return
        if any(abs(x - u) < 1e-9 and abs(y - v) < 1e-9 for u, v, _ in pool): return
        pool.append((x, y, ex))
    lines = list({(b[1] - a[1], a[0] - b[0], (b[1] - a[1]) * a[0] + (a[0] - b[0]) * a[1]) for a, b in itertools.combinations(pts, 2)})
    circles = []
    for a, b, c in itertools.combinations(pts, 3):
        (x1, y1), (x2, y2), (x3, y3) = a, b, c
        d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
        if d == 0: continue
        s1 = x1 * x1 + y1 * y1; s2 = x2 * x2 + y2 * y2; s3 = x3 * x3 + y3 * y3
        ux = (s1 * (y2 - y3) + s2 * (y3 - y1) + s3 * (y1 - y2)) / d
        uy = (s1 * (x3 - x2) + s2 * (x1 - x3) + s3 * (x2 - x1)) / d
        circles.append((ux, uy, (x1 - ux) ** 2 + (y1 - uy) ** 2))
    circles = list(set(circles)); random.shuffle(lines); random.shuffle(circles)
    ops = [("ll", p) for p in itertools.combinations(lines, 2)] + [("lc", (l, c)) for l in lines for c in circles] + [("cc", p) for p in itertools.combinations(circles, 2)]
    random.shuffle(ops)
    for kind, (u, w) in ops:
        if len(pool) >= nmax: break
        if kind == "ll":
            q = exact_line_line(u, w)
            if q: add(float(q[0]), float(q[1]), q)
        elif kind == "lc":
            A, B, C = map(float, u); ux, uy, r2 = map(float, w)
            n = math.hypot(A, B); nx, ny = A / n, B / n; dd = C / n
            t = dd - (nx * ux + ny * uy); h2 = r2 - t * t
            if h2 <= 1e-12: continue
            h = math.sqrt(h2); fx, fy = ux + t * nx, uy + t * ny
            add(fx - h * ny, fy + h * nx); add(fx + h * ny, fy - h * nx)
        else:
            x1, y1, r1s = map(float, u); x2, y2, r2s = map(float, w)
            r1, r2 = math.sqrt(r1s), math.sqrt(r2s); dd = math.hypot(x2 - x1, y2 - y1)
            if dd < 1e-12 or dd > r1 + r2 - 1e-9 or dd < abs(r1 - r2) + 1e-9: continue
            a = (r1 * r1 - r2 * r2 + dd * dd) / (2 * dd); h = math.sqrt(max(r1 * r1 - a * a, 0))
            mx, my = x1 + a * (x2 - x1) / dd, y1 + a * (y2 - y1) / dd
            add(mx - h * (y2 - y1) / dd, my + h * (x2 - x1) / dd); add(mx + h * (y2 - y1) / dd, my - h * (x2 - x1) / dd)
    return pool

def circle_pool(k, rng_seed):
    random.seed(rng_seed)
    pts = [(Fr(0), Fr(0))]; t = 1
    while len(pts) < k + 1:
        x = Fr(1 - t * t, 1 + t * t); y = Fr(2 * t, 1 + t * t)
        for p in [(x, y), (-x, -y), (x, -y), (-x, y)]:
            if p not in pts and len(pts) < k + 1: pts.append(p)
        t += 1
    onc = pts[1:]; pool = [(float(x), float(y), (x, y)) for x, y in pts]
    chords = [(b[1] - a[1], a[0] - b[0], (b[1] - a[1]) * a[0] + (a[0] - b[0]) * a[1]) for a, b in itertools.combinations(onc, 2)]
    random.shuffle(chords)
    for l1, l2 in itertools.combinations(chords, 2):
        q = exact_line_line(l1, l2)
        if q is None or abs(q[0]) > 3 or abs(q[1]) > 3: continue
        if any(abs(float(q[0]) - u) < 1e-9 and abs(float(q[1]) - v) < 1e-9 for u, v, _ in pool): continue
        pool.append((float(q[0]), float(q[1]), q))
        if len(pool) > 70: break
    return pool

def verify_hit(entries):
    if all(e[2] is not None for e in entries):
        r = analyse([e[2] for e in entries])
        return f"EXACT: circles={r['circles']} lines={r['lines']} degenerate={r['degenerate']}"
    P = np.array([(e[0], e[1]) for e in entries]); r = count_circles(P, tol=1e-9)
    return f"float recount (tol 1e-9): {r}"

def sa_search(pool, iters, rng_seed, T0=1.5):
    random.seed(rng_seed); m = len(pool); arr = np.array([(x, y) for x, y, _ in pool])
    def score(S):
        r = count_circles(arr[list(S)]); return 999 if r is None else r[0]
    for _ in range(200):
        S = set(random.sample(range(m), 9)); sc = score(S)
        if sc < 999: break
    else: return (999, [])
    best = (sc, sorted(S))
    for it in range(iters):
        T = T0 * (1 - it / iters) + 0.05
        out = random.choice(sorted(S)); inn = random.randrange(m)
        if inn in S: continue
        S2 = (S - {out}) | {inn}; sc2 = score(S2)
        if sc2 <= sc or random.random() < math.exp((sc - sc2) / T):
            S, sc = S2, sc2
            if sc < best[0]: best = (sc, sorted(S))
    return best

if __name__ == "__main__":
    mode = sys.argv[1]; budget = float(sys.argv[2]); t0 = time.time(); overall = 999; runs = 0; hist = {}
    while time.time() - t0 < budget:
        runs += 1
        pool = seed_pool(random.choice([4, 5, 5, 6]), random.choice([30, 45, 60]), runs) if mode == "closure" else circle_pool(random.choice([8, 12, 16, 20]), runs)
        sc, S = sa_search(pool, 3000, 1000 + runs)
        hist[sc] = hist.get(sc, 0) + 1
        entries = [pool[i] for i in S]
        if sc < overall:
            overall = sc
            print(f"run {runs} pool={len(pool)}: new best circles={sc} pts={[(round(x,5),round(y,5)) for x,y,_ in entries]}", flush=True)
        if sc <= 24:
            print("   hit <= 24 -> verification:", verify_hit(entries), flush=True)
    print(f"{mode}: runs={runs} best={overall} histogram of per-run best counts={dict(sorted(hist.items()))}")
