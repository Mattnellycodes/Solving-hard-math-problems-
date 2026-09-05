"""Self-contained verification that 10 real points carry at most 12 three-point lines (t3(10) <= 12, the
orchard value of Burr–Grünbaum–Sloane 1974), i.e. that no 13 triples on 10 points pairwise sharing <= 1 point
are realisable by 10 distinct real points as exactly-collinear triples.
Counting: the leave graph of a 13-packing has 45 - 39 = 6 edges and odd degree at every point, hence degree
sequence (3,1^9): a claw K_{1,3} plus a perfect matching on the other 6 points.  We enumerate ALL labelled
13-packings by CP-SAT, reduce them to isomorphism classes, and for every class either exhibit an (8_3) or a Fano
sub-configuration (=> unrealisable by facts.py F2/F3), or solve the incidence equations exactly (Groebner basis
over Q with a projective frame; distinctness imposed by saturation).
"""
import itertools, sys, time, json
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2')
from model import ClassStore
import sympy as sp

t0 = time.time()
# Every 13-packing has leave graph = claw K_{1,3} + 3 K_2 (degree sequence (3,1^9): 6 edges, all degrees odd).
# Up to relabelling: centre 0, claw leaves 1,2,3, matching {4,5},{6,7},{8,9}.  Enumerate the decompositions of
# K10 minus this leave into 13 triangles by DFS on the first uncovered edge.
LEAVE = {frozenset((0, 1)), frozenset((0, 2)), frozenset((0, 3)), frozenset((4, 5)), frozenset((6, 7)), frozenset((8, 9))}
EDGES = [frozenset(e) for e in itertools.combinations(range(10), 2) if frozenset(e) not in LEAVE]
assert len(EDGES) == 39
fams = []
def dfs(uncov, chosen):
    if not uncov:
        fams.append(frozenset(chosen)); return
    e = min(uncov, key=lambda f: tuple(sorted(f)))
    a, b = sorted(e)
    for c in range(10):
        if c in e: continue
        T = frozenset((a, b, c))
        if frozenset((a, c)) in uncov and frozenset((b, c)) in uncov:
            dfs(uncov - {e, frozenset((a, c)), frozenset((b, c))}, chosen + [T])
dfs(frozenset(EDGES), [])
print(f"13-triple packings with the fixed leave: {len(fams)}  [{time.time()-t0:.0f}s]", flush=True)
store = ClassStore(10)
for fam in fams:
    store.add(list(fam), [])
print(f"isomorphism classes: {len(store.classes)}  [{time.time()-t0:.0f}s]", flush=True)
def has_sub(fam, r, k):
    """some r-subset carries k triples of fam."""
    return any(sum(1 for S in fam if S <= set(sub)) >= k for sub in itertools.combinations(range(10), r))

def realisable_over_R(lines, charts=None):
    """exact test: coordinates from a projective frame (4 points, no 3 on a line), remaining points constructed
    incrementally (intersection of two known lines / one parameter on a known line / two free parameters),
    Groebner basis of the residual incidence equations saturated by distinctness."""
    lines = [frozenset(l) for l in lines]
    pts = list(range(10))
    frame = next(fr for fr in itertools.combinations(pts, 4) if not any(len(l & set(fr)) >= 3 for l in lines))
    coords = {frame[0]: sp.Matrix([1, 0, 0]), frame[1]: sp.Matrix([0, 1, 0]), frame[2]: sp.Matrix([0, 0, 1]), frame[3]: sp.Matrix([1, 1, 1])}
    params = []; eqs = []; free2 = []
    charts = charts or {}
    def known(l): return [p for p in l if p in coords]
    while len(coords) < 10:
        best = None
        for p in pts:
            if p in coords: continue
            L2 = [l for l in lines if p in l and len(known(l - {p})) >= 2]
            L1 = [l for l in lines if p in l and len(known(l - {p})) == 1]
            sc = (len(L2), len(L1))
            if best is None or sc > best[0]: best = (sc, p, L2, L1)
        sc, p, L2, L1 = best
        if len(L2) >= 2:
            a, b = known(L2[0] - {p})[:2]; c, d = known(L2[1] - {p})[:2]
            coords[p] = (coords[a].cross(coords[b])).cross(coords[c].cross(coords[d]))
        elif len(L2) == 1:
            a, b = known(L2[0] - {p})[:2]; t = sp.Symbol(f't{p}'); params.append(t)
            coords[p] = coords[a] + t * coords[b]       # chart missing p = b (excluded by distinctness)
        else:
            u, v = sp.symbols(f'u{p} v{p}'); params += [u, v]; free2.append(p)
            ch = charts.get(p, 0)                        # the three affine charts (u,v,1), (u,1,v), (1,u,v) cover RP^2
            coords[p] = sp.Matrix([[u, v, 1], [u, 1, v], [1, u, v]][ch])
    for l in lines:
        a, b, c = sorted(l)
        e = sp.expand(sp.Matrix.hstack(coords[a], coords[b], coords[c]).det())
        if e != 0: eqs.append(e)
    # Solve the residual system without saturation (it is tiny), then discard solutions with coincident points.
    G = sp.groebner(eqs, *params, order='lex')
    unit = list(G.exprs) == [1]
    sols = [] if unit else sp.solve(list(G.exprs), params, dict=True)
    # A solution is degenerate if two points coincide or a constructed point is the zero vector (the two lines
    # defining it coincide, i.e. two of the 13 lines are equal).  On a positive-dimensional component we test
    # whether some pair coincides / some point vanishes IDENTICALLY in the remaining parameters.
    good = []
    for s in sols:
        P = {p: sp.simplify(coords[p].subs(s)) for p in pts}
        coinc = [(i, j) for i, j in itertools.combinations(pts, 2) if sp.Matrix.hstack(P[i], P[j]).rank() < 2]
        zero = [p for p in pts if all(sp.simplify(c) == 0 for c in P[p])]
        degenerate = bool(coinc or zero)
        if not degenerate:
            good.append((s, 'zero-dimensional' if len(s) == len(params) else 'POSITIVE-DIMENSIONAL, generic point non-degenerate'))
        else:
            print(f"      solution {s}: degenerate (identically coincident pairs {coinc[:4]}..., zero vectors {zero})", flush=True)
    return unit, params, eqs, G, free2, sols, good

for k, c in enumerate(store.classes):
    fam = c['F']
    mk = has_sub(fam, 8, 8); fano = has_sub(fam, 7, 7)
    line = f"class {k} (copies {c['copies']}): contains (8_3): {mk}, contains Fano: {fano}"
    if not (mk or fano):
        unit, params, eqs, G, free2, sols, good = realisable_over_R(fam)
        line += f"; exact: {len(params)} parameters ({len(free2)} two-parameter points), {len(eqs)} equations, Groebner basis = (1): {unit}; solutions {sols}; solutions with distinct points: {good}"
        ok = unit or not good
        for combo in itertools.product(range(3), repeat=len(free2)):
            if all(c == 0 for c in combo): continue
            unit2, _, _, G2, _, sols2, good2 = realisable_over_R(fam, charts=dict(zip(free2, combo)))
            line += f"\n   charts {dict(zip(free2, combo))}: basis = (1): {unit2}; solutions {sols2}; with distinct points: {good2}"
            ok &= (unit2 or not good2)
        line += f"\n   => {'NOT REALISABLE over C (hence over R)' if ok else 'NOT DECIDED / REALISABLE?'}"
    print(line, flush=True)
    print("   triples:", sorted(sorted(S) for S in fam))
print(f"done [{time.time()-t0:.0f}s]")
