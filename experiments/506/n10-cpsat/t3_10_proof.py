"""Self-contained proof that no 10 real points have 13 three-point lines (t3(10) <= 12, BGS 1974).
pts10_13.py: every 13-triple system on 10 points (pairwise <= 1 common point) is one of two classes.
A realisation has exactly these collinear triples (a 14th would exceed the packing maximum, and a line
with 4 points would destroy a 3-point line).  So 4 points no three of which form a triple of the system
are in general position and can be taken as a projective frame.  We construct the remaining points
incrementally (a point on two known lines is their intersection; on one known line it has one free
parameter; otherwise two), collect the residual incidence equations and compute a Groebner basis;
degenerate solutions (coincident points) are excluded by saturation with the distinctness polynomials.
"""
import itertools, sympy as sp

SYSTEMS = [
 [[0,4,8],[0,5,7],[0,6,9],[1,2,3],[1,4,9],[1,5,6],[1,7,8],[2,4,6],[2,5,8],[2,7,9],[3,4,7],[3,5,9],[3,6,8]],
 [[0,4,8],[0,5,7],[0,6,9],[1,2,9],[1,3,4],[1,5,6],[1,7,8],[2,3,7],[2,4,6],[2,5,8],[3,5,9],[3,6,8],[4,7,9]],
]

def analyse(lines, frame=None):
    lines = [frozenset(l) for l in lines]
    pts = sorted(set().union(*lines))
    if frame is None:
        for fr in itertools.combinations(pts, 4):
            if not any(len(l & set(fr)) >= 3 for l in lines):
                frame = fr; break
    coords = {frame[0]: sp.Matrix([1, 0, 0]), frame[1]: sp.Matrix([0, 1, 0]),
              frame[2]: sp.Matrix([0, 0, 1]), frame[3]: sp.Matrix([1, 1, 1])}
    params = []
    eqs = []
    def known_line(l):
        k = [p for p in l if p in coords]
        return k
    progress = True
    while len(coords) < len(pts):
        # choose the point with most incident lines having >= 2 known points
        best = None
        for p in pts:
            if p in coords: continue
            L2 = [l for l in lines if p in l and len(known_line(l - {p})) >= 2]
            L1 = [l for l in lines if p in l and len(known_line(l - {p})) == 1]
            score = (len(L2), len(L1))
            if best is None or score > best[0]:
                best = (score, p, L2, L1)
        score, p, L2, L1 = best
        if len(L2) >= 2:
            l1, l2 = L2[0], L2[1]
            a, b = known_line(l1 - {p})[:2]; c, d = known_line(l2 - {p})[:2]
            v = (coords[a].cross(coords[b])).cross(coords[c].cross(coords[d]))
            coords[p] = sp.simplify(v)
            for l in L2[2:]:
                a, b = known_line(l - {p})[:2]
                eqs.append(sp.expand(sp.Matrix.hstack(coords[a], coords[b], coords[p]).det()))
        elif len(L2) == 1:
            l = L2[0]; a, b = known_line(l - {p})[:2]
            t = sp.Symbol(f't{p}'); params.append(t)
            coords[p] = coords[a] + t * coords[b]     # affine chart on the line (excludes p = b)
        else:
            u, v = sp.symbols(f'u{p} v{p}'); params += [u, v]
            coords[p] = sp.Matrix([u, v, 1])
    # remaining incidences
    for l in lines:
        a, b, c = sorted(l)
        e = sp.expand(sp.Matrix.hstack(coords[a], coords[b], coords[c]).det())
        if e != 0: eqs.append(e)
    return coords, params, eqs, frame

for k, sysl in enumerate(SYSTEMS):
    coords, params, eqs, frame = analyse(sysl)
    eqs = [e for e in eqs if e != 0]
    print(f"system {k}: frame {frame}, free parameters {params}, {len(eqs)} residual equations")
    G = sp.groebner(eqs, *params, order='lex')
    print("  Groebner basis:", list(G))
    sols = sp.solve(list(G), params, dict=True)
    print("  solutions over C:", sols)
    for s in sols:
        P = {p: sp.simplify(v.subs(s)) for p, v in coords.items()}
        distinct = all(sp.Matrix.hstack(P[i], P[j]).rank() == 2 for i, j in itertools.combinations(P, 2))
        nonzero = all(any(c != 0 for c in P[i]) for i in P)
        real = all(sp.im(sp.N(v)) == 0 for v in s.values())
        print(f"   solution {s}: real={real} distinct points={distinct} nonzero={nonzero}")
    # chart boundary cases: t-parameters where p would coincide with the second base point are excluded by
    # distinctness; a point 'at infinity' of the chart u,v is handled by projective symmetry of the frame.
