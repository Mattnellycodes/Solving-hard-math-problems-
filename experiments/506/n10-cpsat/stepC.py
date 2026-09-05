"""Step C.  By Step A every realisation of the skeleton-3 block structure is Möbius-equivalent to
S = {w^k} (regular pentagon on the unit circle, w = exp(2 pi i/5)) and R = {s*lam*w^k}, s = ±1, lam > 0,
lam != 1.  For an abstract structure (blocks, lines) we (i) find all isomorphisms of its block structure onto
the pentagon block structure (there may be several; each gives a concrete line pattern), (ii) for each
pattern require the point Z (= infinity of the Euclidean picture) to lie on the circle through each line's
three points: with u = Z, v = conj(Z) as independent unknowns,
   det[[u v, u, v, 1], [p pbar, p, pbar, 1], [q qbar, q, qbar, 1], [r rbar, r, rbar, 1]] = 0,
polynomial in (u, v, lam, s, w) with Phi_5(w) = 0, s^2 = 1.  A Groebner basis equal to {1} proves that no
such Z exists for any lam, i.e. the abstract structure is not realisable.  (A line of size 4 = a 4-block
whose circle passes through Z.)  Usage: python3 stepC.py file.json
"""
import sys, json, itertools, math, sympy as sp
import networkx as nx
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat')
from n10model import incidence_graph

u, v, lam, s, w = sp.symbols('u v lam s w')
Phi5 = w**4 + w**3 + w**2 + w + 1
# pentagon labels: S_k = w^k (k = 0..4), R_k = s*lam*w^k (labels 5..9); conj(w) = w^4 (mod Phi5), conj(lam)=lam, conj(s)=s
def coord(k):
    return (w**k, w**(4 * k)) if k < 5 else (s * lam * w**(k - 5), s * lam * w**(4 * (k - 5)))
def circle_condition(p, q, r):
    rows = [[u * v, u, v, 1]]
    for k in (p, q, r):
        z, zb = coord(k); rows.append([z * zb, z, zb, 1])
    return sp.expand(sp.Matrix(rows).det())
def reduce_w(e):
    return sp.rem(sp.expand(e), Phi5, w)
# pentagon block structure (labelled): blocks {a,b,c,d}: a,b in S, c,d in R with chord ab parallel to chord cd:
# chord w^i w^j is parallel to chord w^k w^l  iff  i + j = k + l (mod 5)
def pentagon_blocks():
    F = [frozenset(range(5)), frozenset(range(5, 10))]
    for i, j in itertools.combinations(range(5), 2):
        for k, l in itertools.combinations(range(5), 2):
            if (i + j - k - l) % 5 == 0:
                F.append(frozenset({i, j, 5 + k, 5 + l}))
    return F
PF = pentagon_blocks()
assert len(PF) == 22
GP = incidence_graph(PF, [])
nm = nx.algorithms.isomorphism.categorical_node_match('c', 0)
data = json.load(open(sys.argv[1]))
for idx, rec in enumerate(data):
    F = [frozenset(B) for B in rec['blocks']]; L = [frozenset(S) for S in rec['lines']]
    G = incidence_graph(F, [])
    GM = nx.algorithms.isomorphism.GraphMatcher(G, GP, node_match=nm)
    patterns = set()
    n_iso = 0
    for iso in GM.isomorphisms_iter():
        n_iso += 1
        pm = {p: q[1] for p, q in iso.items() if p[0] == 'p'}
        patterns.add(frozenset(frozenset(pm[x] for x in S_) for S_ in L))
    print(f"structure {idx}: block structure isomorphic to pentagon structure: {n_iso > 0} ({n_iso} isomorphisms, {len(patterns)} distinct line patterns)", flush=True)
    for pat in patterns:
        eqs = [Phi5, s**2 - 1]
        for S_ in pat:
            S_ = sorted(S_)
            for trip in itertools.combinations(S_, 3):      # all triples of a line's points must be concyclic with Z
                eqs.append(reduce_w(circle_condition(*trip)))
        Gb = sp.groebner(eqs, u, v, lam, s, w, order='grevlex')
        unit = (list(Gb) == [1])
        print(f"   pattern {[sorted(S_) for S_ in sorted(pat, key=sorted)]}: Groebner basis = {{1}}: {unit}", flush=True)
        if not unit:
            print("      basis:", list(Gb)[:6], flush=True)
