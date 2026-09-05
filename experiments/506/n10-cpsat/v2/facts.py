"""Self-contained verification of the small combinatorial/geometric facts used by model.py and filters.py
(n10-cpsat/v2 agent).  Run: python3 facts.py
 F1  packing bounds on 9 points (derived structure at a point of a 10-set): <= 2 five-lines, <= 3 four-lines,
     and with k four-lines at most 12, 8, 7, 4 three-lines (k = 0..3), all lines pairwise sharing <= 1 point.
 F2  7 triples on 7 points pairwise sharing <= 1 point = Fano plane (covers all 21 pairs) -> excluded by SG.
 F3  8 triples on 8 points pairwise sharing <= 1 point: unique up to isomorphism (Möbius–Kantor (8_3)), and it has
     no realisation by 8 distinct points of RP^2 (Groebner basis: the parameter satisfies t^2 - t + 1 = 0).
 F4  11 triples on 9 points pairwise sharing <= 1 point: every such system is AG(2,3) minus one line (10080
     labelled systems = 840 x 12) and contains an (8_3) on 8 of its points -> 9 real points carry <= 10
     three-point lines (t3(9) <= 10), so <= 10 four-blocks through any point of a 10-set.
 F5  Miquel: for 8 distinct points of C labelled by the cube vertices, a fixed product of the six face
     cross-ratios (with exponents +-1) is identically 1; hence if five faces are concyclic (real cross-ratio)
     so is the sixth.  The identity is checked symbolically.
"""
import itertools, math, time
from ortools.sat.python import cp_model
import sympy as sp

def packing_max(npts, fixed_sizes, free_size):
    """max number of free_size-subsets pairwise (and with the fixed lines) sharing <= 1 point, given a list of
    fixed line sizes to be placed too (all lines pairwise share <= 1 point)."""
    pts = range(npts)
    m = cp_model.CpModel()
    cand = {}
    for k in set(list(fixed_sizes) + [free_size]):
        for c in itertools.combinations(pts, k):
            cand[frozenset(c)] = m.NewBoolVar('')
    for Q in itertools.combinations(pts, 2):
        m.AddAtMostOne([v for S, v in cand.items() if set(Q) <= S])
    for k in set(fixed_sizes):
        m.Add(sum(v for S, v in cand.items() if len(S) == k) == fixed_sizes.count(k))
    m.Maximize(sum(v for S, v in cand.items() if len(S) == free_size))
    s = cp_model.CpSolver(); s.parameters.num_workers = 2; s.parameters.max_time_in_seconds = 300
    st = s.Solve(m)
    assert st == cp_model.OPTIMAL, s.StatusName(st)
    return int(s.ObjectiveValue())

def all_packings(npts, k, count):
    """all labelled families of `count` k-subsets of range(npts) pairwise sharing <= 1 point."""
    pts = range(npts)
    m = cp_model.CpModel()
    cand = {frozenset(c): m.NewBoolVar('') for c in itertools.combinations(pts, k)}
    for Q in itertools.combinations(pts, 2):
        m.AddAtMostOne([v for S, v in cand.items() if set(Q) <= S])
    m.Add(sum(cand.values()) == count)
    sols = []
    class CB(cp_model.CpSolverSolutionCallback):
        def on_solution_callback(self):
            sols.append(frozenset(S for S, v in cand.items() if self.Value(v)))
    s = cp_model.CpSolver(); s.parameters.enumerate_all_solutions = True; s.parameters.num_workers = 1
    st = s.Solve(m, CB())
    assert st == cp_model.OPTIMAL, s.StatusName(st)
    return sols

def iso_classes(npts, fams):
    """isomorphism classes of families of subsets under S_npts (brute force canonical form; npts <= 9)."""
    reps = {}
    perms = list(itertools.permutations(range(npts)))
    for fam in fams:
        best = None
        for perm in perms:
            img = tuple(sorted(tuple(sorted(perm[i] for i in S)) for S in fam))
            if best is None or img < best:
                best = img
        reps[best] = reps.get(best, 0) + 1
    return reps

if __name__ == '__main__':
    t0 = time.time()
    # ---- F1
    print("F1 packing bounds on 9 points:")
    print("   max five-lines:", packing_max(9, [], 5), " max four-lines:", packing_max(9, [], 4))
    for k in range(4):
        print(f"   with {k} four-lines: max three-lines = {packing_max(9, [4]*k, 3)}")
    # ---- F2
    fam7 = all_packings(7, 3, 7)
    cl7 = iso_classes(7, fam7)
    print(f"F2: 7-triple packings on 7 points: {len(fam7)} labelled, {len(cl7)} class(es); all cover all 21 pairs: "
          f"{all(len({frozenset(Q) for S in fam for Q in itertools.combinations(S, 2)}) == 21 for fam in fam7)}")
    # ---- F3
    fam8 = all_packings(8, 3, 8)
    MK = [frozenset({i % 8, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
    assert frozenset(MK) in fam8
    autMK = sum(1 for perm in itertools.permutations(range(8)) if frozenset(frozenset(perm[i] for i in S) for S in MK) == frozenset(MK))
    print(f"F3: 8-triple packings on 8 points: {len(fam8)} labelled; |Aut(MK)| = {autMK}; 8!/|Aut| = {math.factorial(8)//autMK}"
          f" -> single isomorphism class: {len(fam8) == math.factorial(8)//autMK}  [{time.time()-t0:.0f}s]")
    # realisability of MK over R: frame 0,1,2,5 (no three on a line of MK)
    lines = MK
    assert not any(len(l & {0, 1, 2, 5}) >= 3 for l in lines)
    t = sp.Symbol('t')
    P = {0: sp.Matrix([1, 0, 0]), 1: sp.Matrix([0, 1, 0]), 2: sp.Matrix([0, 0, 1]), 5: sp.Matrix([1, 1, 1])}
    def line(a, b): return P[a].cross(P[b])
    P[3] = line(0, 1).cross(line(2, 5))          # {0,1,3}, {2,3,5}
    # p4 on line(1,2) [{1,2,4}]: one parameter; chart P[4] = P[1] + t*P[2] (misses only p4 = p2, excluded since distinct)
    P[4] = P[1] + t * P[2]
    P[6] = line(5, 0).cross(line(3, 4))          # {5,6,0}, {3,4,6}
    P[7] = line(0, 2).cross(line(4, 5))          # {7,0,2}, {4,5,7}
    eqs = []
    for l in lines:
        a, b, c = sorted(l)
        e = sp.expand(sp.Matrix.hstack(P[a], P[b], P[c]).det())
        if e != 0:
            eqs.append(e)
    G = sp.groebner(eqs, t, order='lex')
    print("   MK residual equations:", eqs, " Groebner basis:", list(G))
    sols = sp.solve(list(G)[0], t)
    print("   solutions t:", sols, " real:", [bool(sp.im(s) == 0) for s in sols])
    # the other chart (p4 = p2 excluded by distinctness): so no real realisation with distinct points.
    # ---- F4
    fam9 = all_packings(9, 3, 11)
    print(f"F4: 11-triple packings on 9 points: {len(fam9)} labelled (expected 10080 = 840 STS(9) x 12 lines)  [{time.time()-t0:.0f}s]")
    # leave graph of each is a triangle -> completes to an STS(9); check directly that adding the uncovered pairs' triple works
    ok_ag = True; ok_mk = True
    for fam in fam9:
        cov = {frozenset(Q) for S in fam for Q in itertools.combinations(S, 2)}
        leave = [Q for Q in itertools.combinations(range(9), 2) if frozenset(Q) not in cov]
        tri = frozenset(itertools.chain.from_iterable(leave))
        if not (len(leave) == 3 and len(tri) == 3):
            ok_ag = False
        # (8_3) inside: some 8-subset carrying 8 of the triples
        has = any(sum(1 for S in fam if S <= set(sub)) == 8 for sub in itertools.combinations(range(9), 8))
        if not has:
            ok_mk = False
    print(f"   every 11-packing = STS(9) minus a line: {ok_ag};  every 11-packing contains an (8_3): {ok_mk}")
    print("   => 9 real points carry at most 10 three-point lines (t3(9) <= 10), given F3.")
    # ---- F5 Miquel identity.  For a face with cyclic vertex order (a,b,c,d) the EDGE cross-ratio
    #   cr = (a-b)(c-d) / ((b-c)(d-a))  ( = (a,c;b,d) )  is real iff a,b,c,d are concyclic/collinear.
    # Orientation rule: faces x=const use (y-edges)/(z-edges), y=const use (z-edges)/(x-edges), z=const use
    # (x-edges)/(y-edges).  Every edge then occurs in exactly two faces with exponents +1 and -1, so the product
    # of the six face cross-ratios is +-1 identically; five real cross-ratios force the sixth to be real.
    verts = list(itertools.product((0, 1), repeat=3))
    z = {v: sp.Symbol('z' + ''.join(map(str, v))) for v in verts}
    def edge_prod(edges):
        return sp.Mul(*[(z[u] - z[v]) for u, v in edges])
    def face_cr(i, val):
        f = [v for v in verts if v[i] == val]
        j, k = [c for c in range(3) if c != i]          # j < k; edges in direction j and direction k
        ej = [(u, v) for u in f for v in f if u[k] == v[k] and u[j] == 0 and v[j] == 1]   # j-edges
        ek = [(u, v) for u in f for v in f if u[j] == v[j] and u[k] == 0 and v[k] == 1]   # k-edges
        # rule: x-face: y/z (j=1,k=2 -> ej/ek); y-face: z/x (j=0,k=2 -> ek/ej); z-face: x/y (j=0,k=1 -> ej/ek)
        num, den = (ej, ek) if i != 1 else (ek, ej)
        return edge_prod(num) / edge_prod(den)
    crs = [face_cr(i, val) for i in range(3) for val in (0, 1)]
    prod = sp.simplify(sp.Mul(*crs))
    print("F5: Miquel: product of the six oriented edge cross-ratios =", prod)
    # each factor is a cross-ratio of its face: check against the definition (p,q;r,s) = (p-r)(q-s)/((p-s)(q-r))
    def cr4(p, q, r, s): return (p - r) * (q - s) / ((p - s) * (q - r))
    ok = True
    for i in range(3):
        for val in (0, 1):
            f = [v for v in verts if v[i] == val]
            j, k = [c for c in range(3) if c != i]
            a = [v for v in f if (v[j], v[k]) == (0, 0)][0]; b = [v for v in f if (v[j], v[k]) == (1, 0)][0]
            c = [v for v in f if (v[j], v[k]) == (1, 1)][0]; d = [v for v in f if (v[j], v[k]) == (0, 1)][0]
            expr = face_cr(i, val)
            cands = [cr4(z[a], z[c], z[b], z[d]), 1 / cr4(z[a], z[c], z[b], z[d])]
            ok &= any(sp.simplify(expr - cnd) == 0 for cnd in cands)
    print("   each factor is a cross-ratio (a,c;b,d) or its inverse of its face's cyclic order (a,b,c,d):", ok)
    print(f"done [{time.time()-t0:.0f}s]")
