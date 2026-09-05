"""Exact non-realisability test for 4-block families containing a 6-block A (own code).

Lemma (pencil involution).  Let A be a circle of the sphere carrying the 6-block, and r, r' two points of P off A.
The circles through r, r' are the plane sections through the line l = rr'; each meets the plane of A in a line
through the point X = l ∩ plane(A), so the pairs {a, b} ⊂ A cut out by circles through r, r' are pairs of the
projective involution of the conic A with centre X (Frégier involution; X is not on A since r, r' are not on A).
Hence if the three 4-blocks {r, r', a, b} through {r, r'} pair up the six points of A (a perfect matching M_rr'),
the three pairs are pairs of ONE involution: with homogeneous coordinates (x_i : y_i) of the six points on the conic
(A = P^1(R)), an involution is  alpha x x' + beta (x y' + x' y) + gamma y y' = 0, so
      det [[x_a x_b, x_a y_b + x_b y_a, y_a y_b] for the three pairs] = 0.
Normalisation: PGL(2,R) is sharply 3-transitive, so t_0 = 0, t_1 = 1, t_2 = infinity; t_3, t_4, t_5 are then finite,
pairwise distinct and different from 0, 1.  A Groebner basis (1) of the determinant conditions saturated by the
distinctness product proves that the 4-block family (whatever the lines) is not realisable.
Usage: python3 analyse_six.py report.json [--all]
"""
import sys, json, itertools, argparse, time
import sympy as sp

ap = argparse.ArgumentParser(); ap.add_argument('file'); ap.add_argument('--all', action='store_true')
args = ap.parse_args()
data = json.load(open(args.file))
x3, x4, x5, w = sp.symbols('x3 x4 x5 w')
seen = set()
for rec in data:
    F = [frozenset(B) for B in rec['blocks']]
    sixes = [B for B in F if len(B) == 6]
    if len(sixes) != 1: continue
    if not args.all and rec['report']['kills']: continue
    key = tuple(sorted(tuple(sorted(B)) for B in F))
    if key in seen:
        print(f"structure {rec['index']}: same 4-block family as an earlier structure"); continue
    seen.add(key)
    A = sorted(sixes[0]); R = [p for p in range(10) if p not in A]
    lab = {A[0]: (0, 1), A[1]: (1, 1), A[2]: (1, 0), A[3]: (x3, 1), A[4]: (x4, 1), A[5]: (x5, 1)}
    conds = []; info = []
    for r, r2 in itertools.combinations(R, 2):
        pairs = [tuple(sorted(B & set(A))) for B in F if {r, r2} <= B and len(B) == 4 and len(B & set(A)) == 2]
        if len(pairs) == 3 and len(set(itertools.chain.from_iterable(pairs))) == 6:
            rows = []
            for a, b in pairs:
                xa, ya = lab[a]; xb, yb = lab[b]
                rows.append([xa * xb, xa * yb + xb * ya, ya * yb])
            conds.append(sp.expand(sp.Matrix(rows).det()))
            info.append(((r, r2), pairs))
    print(f"structure {rec['index']} (count {rec['count']}): 6-block {A}, R = {R}; full matchings: {len(info)}")
    for (rr, pairs) in info:
        print(f"      pair {rr}: matching {pairs}")
    dist = sp.Integer(1)
    for a, b in itertools.combinations(A, 2):
        xa, ya = lab[a]; xb, yb = lab[b]
        dist *= (xa * yb - xb * ya)
    t0 = time.time()
    G = sp.groebner(conds + [sp.expand(dist * w - 1)], w, x3, x4, x5, order='grevlex')
    unit = list(G.exprs) == [1]
    print(f"   {len(conds)} involution conditions; saturated Groebner basis = (1): {unit}  [{time.time()-t0:.1f}s]")
    if not unit:
        print("   basis:", list(G.exprs)[:8])
        sols = sp.solve(conds, [x3, x4, x5], dict=True)
        print("   solutions of the unsaturated system:", sols[:10])
    print("   =>", "4-BLOCK FAMILY NOT REALISABLE (regardless of lines)" if unit else "NOT DECIDED")
