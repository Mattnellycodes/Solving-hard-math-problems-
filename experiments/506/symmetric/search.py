"""Symmetric-orbit search for Erdos #506: unions of D_m orbits with continuous parameters.
Special parameter values are found as roots of coincidence conditions (concyclic quadruples,
collinear triples, quadruples including the centre), one orbit representative per condition.

Usage: python3 search.py [m_min m_max] [--nmax 25] [--out results_1d.json]
"""
import sys, json, math, itertools, time, argparse, os
import numpy as np
import symcore as S

HERE = os.path.dirname(os.path.abspath(__file__))


def gen_batch(m, orbits, params_batch):
    """params_batch: (T, p) -> pts (T, n, 2) (finite points only)."""
    pb = np.asarray(params_batch, float)
    T = pb.shape[0]
    cols = []; k = 0
    for ob in orbits:
        t = ob[0]
        if t == 'A':
            r = pb[:, k]; k += 1
            a = 2 * np.pi * np.arange(m) / m
            cols.append(np.stack([r[:, None] * np.cos(a)[None, :], r[:, None] * np.sin(a)[None, :]], -1))
        elif t == 'B':
            r = pb[:, k]; k += 1
            a = (2 * np.arange(m) + 1) * np.pi / m
            cols.append(np.stack([r[:, None] * np.cos(a)[None, :], r[:, None] * np.sin(a)[None, :]], -1))
        elif t == 'G':
            r, th = pb[:, k], pb[:, k + 1]; k += 2
            base = 2 * np.pi * np.arange(m) / m
            ang = np.stack([th[:, None] + base[None, :], -th[:, None] + base[None, :]], -1).reshape(T, 2 * m)
            cols.append(np.stack([r[:, None] * np.cos(ang), r[:, None] * np.sin(ang)], -1))
        elif t == 'R':
            r, th = pb[:, k], pb[:, k + 1]; k += 2
            base = 2 * np.pi * np.arange(m) / m
            ang = th[:, None] + base[None, :]
            cols.append(np.stack([r[:, None] * np.cos(ang), r[:, None] * np.sin(ang)], -1))
        elif t == 'C':
            cols.append(np.zeros((T, 1, 2)))
    return np.concatenate(cols, axis=1)


class Family:
    def __init__(self, m, orbits):
        self.m = m; self.orbits = orbits
        self.has_inf = any(o[0] == 'I' for o in orbits)
        self.has_c = any(o[0] == 'C' for o in orbits)
        self.perms, self.labels = S.group_perms(m, orbits)
        self.n = len(self.labels)
        self.N = self.n + (1 if self.has_inf else 0)
        self.p = S.nparams(orbits)
        # parameter kinds: 'r' radius, 'a' angle
        self.kinds = []
        for o in orbits:
            if o[0] in 'AB': self.kinds.append('r')
            elif o[0] == 'G': self.kinds += ['r', 'a']
            elif o[0] == 'R': self.kinds += ['r', 'f']     # 'f' = full-range rotation angle
        # extended index set: finite points + virtual centre (if not present)
        n_ext = self.n + (0 if self.has_c else 1)
        self.n_ext = n_ext
        perms_ext = [np.concatenate([p, np.arange(self.n, n_ext)]) for p in self.perms]
        self.quads = np.array(S.orbit_reps(itertools.combinations(range(n_ext), 4), perms_ext)).reshape(-1, 4)
        self.trips = np.array(S.orbit_reps(itertools.combinations(range(n_ext), 3), perms_ext)).reshape(-1, 3)
        self.K = len(self.quads) + len(self.trips)

    def name(self):
        return f"D{self.m} " + " ".join(o[0] for o in self.orbits)

    def points(self, params):
        return gen_batch(self.m, self.orbits, np.asarray(params, float)[None, :])[0]

    def ext_batch(self, pb):
        P = gen_batch(self.m, self.orbits, pb)
        if not self.has_c:
            P = np.concatenate([P, np.zeros((P.shape[0], 1, 2))], axis=1)
        return P

    def conds_batch(self, pb):
        """(T, K) array of condition values (normalised)."""
        P = self.ext_batch(pb)
        T = P.shape[0]
        Q = P[:, self.quads]          # (T, kq, 4, 2)
        x = Q[..., 0]; y = Q[..., 1]
        M = np.stack([x, y, x * x + y * y, np.ones_like(x)], -1)
        dq = np.linalg.det(M)
        # normalise by product of scales to make tolerances meaningful
        sc = np.maximum(np.abs(Q).max(axis=(2, 3)), 1e-9)
        dq = dq / sc ** 4
        Tq = P[:, self.trips]
        x = Tq[..., 0]; y = Tq[..., 1]
        dt = (x[..., 1] - x[..., 0]) * (y[..., 2] - y[..., 0]) - (x[..., 2] - x[..., 0]) * (y[..., 1] - y[..., 0])
        sc = np.maximum(np.abs(Tq).max(axis=(2, 3)), 1e-9)
        dt = dt / sc ** 2
        return np.concatenate([dq, dt], axis=1)

    def cond_single(self, pb, cols):
        """values of condition cols[t] at params pb[t]."""
        P = self.ext_batch(pb)
        out = np.empty(len(cols))
        nq = len(self.quads)
        cols = np.asarray(cols)
        isq = cols < nq
        if isq.any():
            idx = np.where(isq)[0]
            Q = P[idx[:, None], self.quads[cols[idx]]]      # (k,4,2)
            x = Q[..., 0]; y = Q[..., 1]
            M = np.stack([x, y, x * x + y * y, np.ones_like(x)], -1)
            out[idx] = np.linalg.det(M) / np.maximum(np.abs(Q).max(axis=(1, 2)), 1e-9) ** 4
        if (~isq).any():
            idx = np.where(~isq)[0]
            Tq = P[idx[:, None], self.trips[cols[idx] - nq]]
            x = Tq[..., 0]; y = Tq[..., 1]
            out[idx] = ((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])) / np.maximum(np.abs(Tq).max(axis=(1, 2)), 1e-9) ** 2
        return out

    def grid(self, kind, G):
        if kind == 'r':
            return np.exp(np.linspace(math.log(0.08), math.log(12.5), G))
        if kind == 'f':
            return np.linspace(0.003, 2 * math.pi / self.m - 0.003, G)
        return np.linspace(0.003, math.pi / self.m - 0.003, G)

    def roots_along(self, base, i, G=700, tol=1e-12):
        """All parameter values t of parameter i (others fixed at base) where some condition vanishes."""
        base = np.asarray(base, float)
        grid = self.grid(self.kinds[i], G)
        pb = np.repeat(base[None, :], G, 0); pb[:, i] = grid
        V = self.conds_batch(pb)                        # (G,K)
        sg = np.sign(V)
        ch = np.where((sg[:-1] * sg[1:]) < 0)          # sign changes
        lo = grid[ch[0]].copy(); hi = grid[ch[0] + 1].copy(); cols = ch[1]
        flo = V[ch[0], ch[1]].copy()
        # exclude conditions that vanish identically (|V| tiny everywhere)
        ident = np.abs(V).max(0) < 1e-9
        keep = ~ident[cols]
        lo, hi, cols, flo = lo[keep], hi[keep], cols[keep], flo[keep]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            pb = np.repeat(base[None, :], len(mid), 0); pb[:, i] = mid
            fm = self.cond_single(pb, cols)
            left = (np.sign(fm) == np.sign(flo)) | (fm == 0)
            lo = np.where(left, mid, lo); flo = np.where(left, fm, flo)
            hi = np.where(left, hi, mid)
            if np.all(hi - lo < tol): break
        roots = list(0.5 * (lo + hi))
        # also near-tangential zeros: local minima of |V| that are tiny
        A = np.abs(V)
        for c in range(V.shape[1]):
            if ident[c]: continue
            col = A[:, c]
            loc = np.where((col[1:-1] < col[:-2]) & (col[1:-1] < col[2:]) & (col[1:-1] < 1e-4))[0] + 1
            for g in loc:
                if np.sign(V[g - 1, c]) != np.sign(V[g + 1, c]): continue   # handled by bisection
                # refine by golden section on |f|
                a, b = grid[g - 1], grid[g + 1]
                for _ in range(80):
                    m1 = a + 0.382 * (b - a); m2 = a + 0.618 * (b - a)
                    pb = np.repeat(base[None, :], 2, 0); pb[0, i] = m1; pb[1, i] = m2
                    f = np.abs(self.cond_single(pb, [c, c]))
                    if f[0] < f[1]: b = m2
                    else: a = m1
                t = 0.5 * (a + b)
                pb = np.repeat(base[None, :], 1, 0); pb[0, i] = t
                if abs(self.cond_single(pb, [c])[0]) < 1e-9: roots.append(t)
        roots = sorted(set(round(float(t), 12) for t in roots))
        return roots


def nice_values(kind, m):
    if kind == 'f':
        return [math.pi / (2 * m), math.pi / (3 * m), 2 * math.pi / (3 * m), math.pi / (4 * m), math.pi / (6 * m), 4 * math.pi / (3 * m)]
    if kind == 'r':
        return [2.0, 3.0, 0.5, 1 / math.cos(math.pi / m) if m >= 3 else 1.5, (1 + math.sqrt(5)) / 2, math.sqrt(2), math.sqrt(3), 1.3, 2.5, 4.0]
    return [math.pi / (2 * m), math.pi / (3 * m), math.pi / (4 * m), 2 * math.pi / (3 * m), math.pi / (6 * m), math.pi / (5 * m)]


def canonical_orbits(orbits):
    key = {'A': 0, 'B': 1, 'G': 2, 'R': 2.5, 'C': 3, 'I': 4}
    obs = sorted(orbits, key=lambda o: key[o[0]])
    if not any(o[0] == 'A' for o in obs):
        obs = [('A',) if o[0] == 'B' else o for o in obs]
        obs = sorted(obs, key=lambda o: key[o[0]])
    return tuple(obs)


def enumerate_families(m_min, m_max, nmin=5, nmax=25, max_orbits=3, types='ABG', require=None):
    fams = set()
    for m in range(m_min, m_max + 1):
        size = {'A': m, 'B': m, 'G': 2 * m, 'R': m}
        for k in range(1, max_orbits + 1):
            for combo in itertools.combinations_with_replacement(types, k):
                if require and require not in combo: continue
                for extra in [(), ('C',), ('I',), ('C', 'I')]:
                    orbs = [(t,) for t in combo] + [(e,) for e in extra]
                    n = sum(size[t] for t in combo) + len(extra)
                    if nmin <= n <= nmax:
                        fams.add((m, canonical_orbits(orbs)))
    return sorted(fams, key=lambda f: (sum({'A': f[0], 'B': f[0], 'G': 2 * f[0], 'R': f[0], 'C': 1, 'I': 1}[o[0]] for o in f[1]), f[0]))


class Results:
    def __init__(self, path):
        self.path = path
        self.best = {}      # n -> list of entries (unique structures)
        self.seen = set()
        if os.path.exists(path):
            d = json.load(open(path))
            self.best = {int(k): v for k, v in d.items()}
            for n, L in self.best.items():
                for e in L: self.seen.add((n, e['sig']))

    def threshold(self, n):
        return min([S.formula(n)] + [e['best'] for e in self.best.get(n, [])])

    def add(self, fam, params, res, note=''):
        n = fam.N
        thr = S.formula(n)
        cur = self.best.get(n, [])
        curbest = min([e['best'] for e in cur], default=10 ** 9)
        if res['best'] > thr or res['best'] > curbest: return False
        sig = f"{res['best']}|{res['nblocks']}|{res['deg']}|{res['lines']}|" + ','.join(map(str, res['sizes']))
        if (n, sig) in self.seen: return False
        self.seen.add((n, sig))
        entry = {'n': n, 'best': res['best'], 'formula': thr, 'family': fam.name(), 'm': fam.m,
                 'orbits': [list(o) for o in fam.orbits], 'params': [float(x) for x in params],
                 'nblocks': res['nblocks'], 'deg': res['deg'], 'lines': res['lines'], 'sizes': res['sizes'],
                 'O': (None if res['O'] in ('inf', None) else [float(res['O'][0]), float(res['O'][1])]),
                 'O_is_inf': res['O'] == 'inf', 'O_blocks': res['O_blocks'], 'sig': sig, 'note': note}
        if res['best'] < curbest:
            self.best[n] = [entry]
        else:
            self.best[n] = cur + [entry]
        self.best[n] = sorted(self.best[n], key=lambda e: e['best'])[:12]
        self.save()
        return True

    def save(self):
        json.dump({str(k): v for k, v in sorted(self.best.items())}, open(self.path, 'w'), indent=1)


def evaluate_params(fam, params):
    pts = fam.points(params)
    return S.evaluate_fast(pts, fam.has_inf)


def search_family(fam, results, G=700, time_budget=120.0, log=print, depth=2):
    t0 = time.time()
    p = fam.p
    # fix the first radius parameter at 1 (scale)
    fixed = None
    for i, k in enumerate(fam.kinds):
        if k == 'r': fixed = i; break
    free = [i for i in range(p) if i != fixed]
    bases = []
    if not free:
        bases.append(np.array([1.0] * p))
    else:
        # base values for the free params: nice values (product), capped
        lists = [nice_values(fam.kinds[i], fam.m) for i in free]
        prods = list(itertools.product(*[range(len(L)) for L in lists]))
        rng = np.random.default_rng(12345)
        if len(prods) > 24: prods = [prods[j] for j in rng.choice(len(prods), 24, replace=False)]
        for pr in prods:
            b = np.ones(p)
            for i, j in zip(free, pr): b[i] = lists[i and free.index(i) or free.index(i)][j] if False else lists[free.index(i)][j]
            bases.append(b)
    n_eval = 0; n_roots = 0
    evaluated = set()
    def try_params(params, note):
        nonlocal n_eval
        key = tuple(round(float(x), 9) for x in params)
        if key in evaluated: return None
        evaluated.add(key)
        r = evaluate_params(fam, params)
        n_eval += 1
        if r is None: return None
        if results.add(fam, params, r, note):
            log(f"  NEW n={fam.N} best={r['best']} (formula {S.formula(fam.N)}) {fam.name()} params={np.round(params,6).tolist()} nb={r['nblocks']} deg={r['deg']} lines={r['lines']} {note}")
        return r
    for b in bases:
        try_params(b, 'base')
    # level-1: 1D root search along each free param from each base; level-2: from each root along other params
    queue = [(b, 0) for b in bases]
    while queue and time.time() - t0 < time_budget:
        base, lvl = queue.pop(0)
        for i in free:
            if time.time() - t0 > time_budget: break
            roots = fam.roots_along(base, i, G=G)
            n_roots += len(roots)
            for t in roots:
                params = base.copy(); params[i] = t
                r = try_params(params, f'root p{i} lvl{lvl}')
                if r is not None and lvl + 1 < depth and len(free) > 1:
                    # continue coordinate search from promising roots only
                    if r['best'] <= S.formula(fam.N) + max(2, fam.N // 2):
                        queue.append((params, lvl + 1))
    log(f"{fam.name()} n={fam.N} p={p} K={fam.K}: {n_roots} roots, {n_eval} evals, {time.time()-t0:.1f}s")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('m_min', type=int, nargs='?', default=2)
    ap.add_argument('m_max', type=int, nargs='?', default=12)
    ap.add_argument('--nmin', type=int, default=5)
    ap.add_argument('--nmax', type=int, default=25)
    ap.add_argument('--out', default=os.path.join(HERE, 'results_1d.json'))
    ap.add_argument('--budget', type=float, default=120.0)
    ap.add_argument('--G', type=int, default=700)
    ap.add_argument('--depth', type=int, default=2)
    ap.add_argument('--max_orbits', type=int, default=3)
    ap.add_argument('--types', default='ABG')
    ap.add_argument('--require', default=None)
    args = ap.parse_args()
    results = Results(args.out)
    fams = enumerate_families(args.m_min, args.m_max, args.nmin, args.nmax, args.max_orbits, args.types, args.require)
    print(f"{len(fams)} families", flush=True)
    for m, orbs in fams:
        fam = Family(m, list(orbs))
        search_family(fam, results, G=args.G, time_budget=args.budget, log=lambda s: print(s, flush=True), depth=args.depth)
    print("SUMMARY")
    for n in sorted(results.best):
        e = results.best[n][0]
        print(f"n={n:2d} formula={S.formula(n):3d} best={e['best']:3d} {e['family']} params={np.round(e['params'],5).tolist()} ({len(results.best[n])} structures)")
