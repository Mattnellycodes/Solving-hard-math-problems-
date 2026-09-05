#!/usr/bin/env python3
"""involution_lemma.py -- exact algebraic verification of the "pencil involution" obstruction for the two
surviving 6-block structures of the n = 10 enumeration.

Setting.  B = {0,..,5} is a 6-block (a circle C on the sphere), a, b are two points off C, and the three blocks
{a,b,x,y} through {a,b} pair the six points of C by the perfect matching M_ab.  The pencil of circles
through a and b cuts C in the pairs of a projective involution of C = P^1(R) (its chords pass through the
radical centre; Fregier involution), so the three pairs of M_ab are three pairs of ONE real projective
involution.  Three pairs (x,x'), (y,y'), (z,z') lie in a common involution  A t t' + B (t + t') + C = 0  iff
        det [[x x', x + x', 1], [y y', y + y', 1], [z z', z + z', 1]] = 0.
Hence a realisation gives six distinct points of P^1(R) satisfying this determinant condition for every
matching that occurs.  We fix three of the points to 0, 1, -1 (PGL(2,R) is sharply 3-transitive) and treat
the case that one of the remaining points is at infinity separately (homogeneous form).  A Groebner basis
of the resulting ideal saturated by the distinctness conditions decides whether COMPLEX solutions exist;
if the saturated ideal is (1) there is no complex solution, a fortiori no real one, and the structure is
not realisable.  (The group-theoretic argument in the report shows more: at most 4 matchings of six
concyclic points can be pencil-induced; here we only need the specific matchings of the two candidates.)
"""
import itertools, sys, json
import sympy as sp


def inv_det(p, q, r, s, t, u):
    """determinant condition for the pairs (p,q), (r,s), (t,u) to be in involution (affine coordinates)."""
    return sp.Matrix([[p * q, p + q, 1], [r * s, r + s, 1], [t * u, t + u, 1]]).det()


def inv_det_h(P, Q, R, S, T, U):
    """homogeneous version: points as (u, v) pairs."""
    def row(A, B):
        return [A[0] * B[0], A[0] * B[1] + A[1] * B[0], A[1] * B[1]]
    return sp.Matrix([row(P, Q), row(R, S), row(T, U)]).det()


def matchings_of(blocks, B, outside):
    """for each outside pair, the perfect matching of B induced by the blocks {a,b,x,y}."""
    res = {}
    for a, b in itertools.combinations(outside, 2):
        M = []
        for bl in blocks:
            if a in bl and b in bl and len(bl) == 4 and len(set(bl) & set(B)) == 2:
                M.append(tuple(sorted(set(bl) - {a, b})))
        assert len(M) == 3 and set(sum(M, ())) == set(B), (a, b, M)
        res[(a, b)] = tuple(sorted(M))
    return res


def check_structure(blocks, verbose=True):
    B = [bl for bl in blocks if len(bl) == 6][0]
    outside = [p for p in range(10) if p not in B]
    ms = matchings_of(blocks, B, outside)
    distinct = sorted(set(ms.values()))
    if verbose:
        for k, v in ms.items():
            print(f"   outside pair {k}: matching {v}")
        print(f"   {len(distinct)} distinct matchings")
    # variables: points 0,1,2 fixed to 0, 1, -1; points 3,4,5 free (affine), plus the cases where one of them is inf
    x3, x4, x5 = sp.symbols('x3 x4 x5')
    results = {}
    for inf_pt in [None, 3, 4, 5]:
        pts = {0: (sp.Integer(0), sp.Integer(1)), 1: (sp.Integer(1), sp.Integer(1)), 2: (sp.Integer(-1), sp.Integer(1)),
               3: (x3, sp.Integer(1)), 4: (x4, sp.Integer(1)), 5: (x5, sp.Integer(1))}
        free = [x3, x4, x5]
        if inf_pt is not None:
            pts[inf_pt] = (sp.Integer(1), sp.Integer(0))
            free = [v for v in free if v != {3: x3, 4: x4, 5: x5}[inf_pt]]
        eqs = []
        for M in distinct:
            (p, q), (r, s), (t, u) = M
            eqs.append(sp.expand(inv_det_h(pts[p], pts[q], pts[r], pts[s], pts[t], pts[u])))
        # distinctness: u_i v_j - u_j v_i != 0 for all pairs -> saturate with a Rabinowitsch variable
        z = sp.symbols('z')
        prod = sp.Integer(1)
        for i, j in itertools.combinations(range(6), 2):
            d = pts[i][0] * pts[j][1] - pts[j][0] * pts[i][1]
            if d.free_symbols:
                prod *= d
        gens = free + [z]
        G = sp.groebner(eqs + [sp.expand(prod * z - 1)], *gens, order='grevlex')
        trivial = (list(G.exprs) == [1])
        if trivial:
            results[inf_pt] = True
            if verbose:
                print(f"   case point-at-infinity={inf_pt}: equations={len(eqs)}, saturated Groebner basis = (1)"
                      f"  -> no complex solution")
        else:
            # zero-dimensional or not: solve exactly and look for real solutions
            sols = sp.solve(list(G.exprs), gens, dict=True)
            real = [s for s in sols if all(sp.im(sp.nsimplify(v)) == 0 and v.is_real for v in s.values())]
            results[inf_pt] = (len(real) == 0)
            if verbose:
                print(f"   case point-at-infinity={inf_pt}: equations={len(eqs)}, saturated Groebner basis = {G.exprs}")
                print(f"      exact solutions: {sols}")
                print(f"      real solutions: {len(real)}  -> {'no real solution' if not real else 'REAL SOLUTION EXISTS'}")
    return all(results.values()), results


if __name__ == "__main__":
    data = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'n10_table_B.json'))
    allok = True
    for i, rec in enumerate(data['results']):
        if 6 not in rec['sizes']:
            continue
        print(f"structure {i}: blocks {rec['blocks']}")
        ok, res = check_structure([tuple(b) for b in rec['blocks']])
        print(f"   => {'NOT REALISABLE over R (no real solution with distinct points)' if ok else 'inconclusive'}")
        allok &= ok
    print("ALL 6-BLOCK CANDIDATES EXCLUDED" if allok else "SOME CASE INCONCLUSIVE")
