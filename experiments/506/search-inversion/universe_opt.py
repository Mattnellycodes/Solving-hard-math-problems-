"""Universe optimisation for Erdos #506.

A *universe* U is a rich Moebius configuration (points on the sphere).  For an n-subset S of U and a
sphere point O not in S, the number of circles of the stereographic projection of S from O equals

    #{ blocks b of U : |b ∩ S| >= 3 and O ∉ b }

because every three points of S span a block of U.  So, once the incidence structure of U is known,
the best planar picture obtainable from ANY n-subset of U and ANY inversion centre that is a point of
U (or an extra high-degree point of the sphere) is an exact combinatorial optimisation problem, which
we solve with CP-SAT for every candidate centre O.

Result: for each universe and n, the minimum count, the subset S and the centre O; every optimum can
then be certified exactly with exact_verify.py (the universes carry exact coordinates when possible).
"""
import numpy as np, itertools, json, time, sys, math
from collections import defaultdict
from ortools.sat.python import cp_model
import sphere as SP


class Universe:
    def __init__(self, name, V, labels=None, exact=None, extra_centres=True, min_extra_deg=3):
        """V: N x 3 unit vectors.  labels: names.  exact: optional dict with exact planar coordinates
        ('planar': list of sympy pairs or 'inf', 'radicands': tuple) for certification."""
        self.name = name
        self.V = SP.normalise_rows(np.asarray(V, float))
        self.N = len(self.V)
        self.labels = labels or [str(i) for i in range(self.N)]
        self.exact = exact
        self.blocks = SP.find_blocks(self.V)
        self.Nn, self.D = SP.planes(self.V, self.blocks)
        self.blocks_of_point = defaultdict(list)
        for t, b in enumerate(self.blocks):
            for p in b:
                self.blocks_of_point[p].append(t)
        self.pdeg = np.array([len(self.blocks_of_point[i]) for i in range(self.N)])
        # extra centres: sphere points on >= min_extra_deg blocks that are not universe points
        self.extras = np.zeros((0, 3)); self.extra_deg = np.zeros(0, int); self.extra_blocks = []
        if extra_centres and len(self.blocks) >= 2:
            C = SP.candidate_points(self.Nn, self.D, max_pairs=4_000_000)
            Uc, mult = SP.merge_points(C, min_mult=2)
            if len(Uc):
                dist = np.min(np.linalg.norm(Uc[:, None, :] - self.V[None, :, :], axis=2), axis=1)
                Uc = Uc[dist > 1e-6]
                deg = SP.degrees(Uc, self.Nn, self.D)
                keep = deg >= min_extra_deg
                self.extras = Uc[keep]; self.extra_deg = deg[keep]
                self.extra_blocks = [SP.blocks_through(x, self.Nn, self.D) for x in self.extras]

    def size_profile(self):
        s = defaultdict(int)
        for b in self.blocks:
            s[len(b)] += 1
        return dict(sorted(s.items()))

    def info(self):
        return dict(name=self.name, N=self.N, nblocks=len(self.blocks), sizes=self.size_profile(),
                    max_pdeg=int(self.pdeg.max()) if self.N else 0,
                    n_extras=len(self.extras),
                    max_extra_deg=int(self.extra_deg.max()) if len(self.extra_deg) else 0)

    def point_classes(self):
        """Crude equivalence classes of points (2 rounds of colour refinement on the incidence graph)."""
        col = [tuple(sorted(len(self.blocks[t]) for t in self.blocks_of_point[i])) for i in range(self.N)]
        for _ in range(2):
            new = []
            for i in range(self.N):
                nb = []
                for t in self.blocks_of_point[i]:
                    nb.append((len(self.blocks[t]), tuple(sorted(col[j] for j in self.blocks[t] if j != i))))
                new.append((col[i], tuple(sorted(nb))))
            # compress
            keys = {k: idx for idx, k in enumerate(sorted(set(new)))}
            col = [keys[k] for k in new]
        classes = defaultdict(list)
        for i, c in enumerate(col):
            classes[c].append(i)
        return list(classes.values())

    # ------------------------------------------------------------------ evaluation of a subset
    def count(self, S, O_index=None, O_extra=None):
        """Exact combinatorial count for subset S (list of indices) and centre (universe point index
        or extra index)."""
        S = set(S)
        if O_index is not None:
            through = set(self.blocks_of_point[O_index])
        else:
            through = set(self.extra_blocks[O_extra])
        c = 0
        for t, b in enumerate(self.blocks):
            if t in through:
                continue
            if len(b & S) >= 3:
                c += 1
        return c

    def blocks_of_subset(self, S):
        S = set(S)
        return [b & S for b in self.blocks if len(b & S) >= 3]


def optimise(U, n, centre, time_limit=60.0, workers=2, forbid=None, log=False):
    """centre = ('point', i) or ('extra', j).  Minimise the number of blocks of S avoiding the centre
    over n-subsets S of U \\ {centre}.  Returns (status, count, S) ."""
    kind, idx = centre
    if kind == 'point':
        through = set(U.blocks_of_point[idx]); banned = {idx}
    else:
        through = set(U.extra_blocks[idx]); banned = set()
    m = cp_model.CpModel()
    x = {p: m.NewBoolVar(f"x{p}") for p in range(U.N) if p not in banned}
    m.Add(sum(x.values()) == n)
    y = {}
    for t, b in enumerate(U.blocks):
        pts = [p for p in b if p in x]
        if len(pts) < 3:
            continue
        # non-degeneracy: not all of S on one block
        if len(pts) >= n:
            m.Add(sum(x[p] for p in pts) <= n - 1)
        if t in through:
            continue
        y[t] = m.NewBoolVar(f"y{t}")
        m.Add(sum(x[p] for p in pts) <= 2 + (len(pts) - 2) * y[t])
    if forbid:
        for S0 in forbid:
            m.Add(sum(x[p] for p in S0 if p in x) <= len(S0) - 1)
    m.Minimize(sum(y.values()))
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_workers = workers
    solver.parameters.log_search_progress = log
    st = solver.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        S = sorted(p for p in x if solver.Value(x[p]))
        return ('OPTIMAL' if st == cp_model.OPTIMAL else 'FEASIBLE', int(round(solver.ObjectiveValue())), S,
                int(round(solver.BestObjectiveBound())))
    return ('INFEASIBLE' if st == cp_model.INFEASIBLE else 'UNKNOWN', None, None, None)


def run_universe(U, ns, time_limit=60.0, workers=2, use_classes=True, out=None, verbose=True,
                 max_centres=None, include_extras=True):
    """For each n, try every centre class representative (universe points) and every extra centre."""
    results = {}
    if use_classes:
        reps = [cl[0] for cl in U.point_classes()]
    else:
        reps = list(range(U.N))
    # order representatives by degree (high degree first)
    reps.sort(key=lambda i: -U.pdeg[i])
    if max_centres:
        reps = reps[:max_centres]
    centres = [('point', i) for i in reps]
    if include_extras:
        order = np.argsort(-U.extra_deg)
        centres += [('extra', int(j)) for j in order[:max_centres or len(order)]]
    for n in ns:
        if n + 1 > U.N:
            continue
        best = None
        t0 = time.time()
        for c in centres:
            st, cnt, S, bound = optimise(U, n, c, time_limit=time_limit, workers=workers)
            if cnt is None:
                continue
            if best is None or cnt < best['count'] or (cnt == best['count'] and st == 'OPTIMAL' and best['status'] != 'OPTIMAL'):
                best = dict(universe=U.name, n=n, count=cnt, status=st, bound=bound, centre=c,
                            centre_label=(U.labels[c[1]] if c[0] == 'point' else f"extra{c[1]}(deg {U.extra_deg[c[1]]})"),
                            S=S, S_labels=[U.labels[p] for p in S])
        if best is not None:
            best['formula'] = SP.formula(n)
            best['seconds'] = round(time.time() - t0, 1)
            # cross-check with direct evaluation
            chk = U.count(best['S'], O_index=best['centre'][1] if best['centre'][0] == 'point' else None,
                          O_extra=best['centre'][1] if best['centre'][0] == 'extra' else None)
            best['check'] = chk
            results[n] = best
            if verbose:
                flag = "  <<< BELOW FORMULA" if best['count'] < best['formula'] else ("  (tie)" if best['count'] == best['formula'] else "")
                print(f"[{U.name}] n={n}: best={best['count']} (f={best['formula']}) status={best['status']} "
                      f"centre={best['centre_label']} S={best['S_labels']} check={chk} t={best['seconds']}s{flag}", flush=True)
            if out:
                with open(out, 'a') as f:
                    f.write(json.dumps({k: (v if not isinstance(v, tuple) else list(v)) for k, v in best.items()}) + "\n")
    return results
