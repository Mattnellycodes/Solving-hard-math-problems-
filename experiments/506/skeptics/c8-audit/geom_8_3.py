"""Independent audit of Lemma 2.4: (a) every family of 8 triples on 8 points pairwise sharing <= 1 point
is isomorphic to the Moebius-Kantor configuration MK = {i,i+1,i+3 mod 8};  (b) MK has no realisation by
8 distinct points of RP^2 (hence none in R^2).  Own derivation: frame {0,1,2,5} (checked to contain no
MK line), p3 forced as the meet of two frame lines, p6/p7/p4 parametrised on their frame lines (the
excluded parameter value is a frame point, impossible for distinct points), then solve exactly."""
import itertools, sympy as sp
n = 8
tri = [frozenset(c) for c in itertools.combinations(range(n), 3)]
perms = list(itertools.permutations(range(n)))
def canon(F):
    return min(tuple(sorted(tuple(sorted(s[i] for i in T)) for T in F)) for s in perms)
# (a) enumerate: any point in >= 4 triples needs >= 8 other points -> impossible; so degrees <= 3 and
# 8 triples * 3 = 24 = 8 * 3 forces all degrees = 3.  We enumerate without using that (only <=1 rule).
classes = set(); count = 0
def rec(F, start):
    global count
    if len(F) == 8:
        count += 1; classes.add(canon(F)); return
    for j in range(start, len(tri)):
        T = tri[j]
        if all(len(T & U) <= 1 for U in F):
            rec(F + [T], j + 1)
# symmetry: the lexicographically first triple {0,1,2} may be assumed present (relabel)
rec([tri[0]], 1)
MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
print(f"labelled 8-triple systems containing {{0,1,2}}: {count}; isomorphism classes: {len(classes)}; "
      f"MK is the class: {canon(MK) in classes}")
# (b) realisation
frame = (0, 1, 2, 5)
assert not any(set(c) <= L for c in itertools.combinations(frame, 3) for L in MK), "frame has a collinear triple"
a, b, c = sp.symbols('a b c')
p = {0: sp.Matrix([1, 0, 0]), 1: sp.Matrix([0, 1, 0]), 2: sp.Matrix([0, 0, 1]), 5: sp.Matrix([1, 1, 1])}
def on_lines(q):  # MK lines through q that contain two frame points
    return [L for L in MK if q in L and len(L & set(frame)) == 2]
# p3: lines {0,1,3} and {2,3,5} -> forced
L1, L2 = on_lines(3)
assert len(on_lines(3)) == 2
def line_through(u, v): return u.cross(v)
def meet(l1, l2): return l1.cross(l2)
p[3] = meet(line_through(p[0], p[1]), line_through(p[2], p[5]))
print("p3 forced =", list(p[3]))
# p6 on line(p5,p0) = s*p5 + t*p0 ; t-free chart: p6 = p5 + a*p0? -> must exclude the frame point itself.
# general point on line(u,v): s*u + t*v.  s = 0 gives v (frame point, excluded), so WLOG s = 1: u + t v.
for q, param in ((6, a), (7, b), (4, c)):
    Ls = on_lines(q); assert len(Ls) == 1, (q, Ls)
    u, v = sorted(Ls[0] & set(frame))
    p[q] = p[u] + param * p[v]
    print(f"p{q} = p{u} + {param}*p{v} =", list(p[q]), f"(excluded value = p{u} itself)")
eqs = [sp.expand(sp.Matrix.hstack(*[p[i] for i in sorted(L)]).det()) for L in MK]
print("incidence equations:", [e for e in eqs if e != 0])
sols = sp.solve([e for e in eqs if e != 0], [a, b, c], dict=True)
print("solutions:", sols)
G = sp.groebner([e for e in eqs if e != 0], a, b, c, order='lex')
uni = [g for g in G if g.free_symbols <= {c}]
print("eliminant in c:", uni, " real roots:", [sp.real_roots(sp.Poly(g, c)) for g in uni])
real = [s for s in sols if all(sp.simplify(sp.im(v)) == 0 for v in s.values())]
print("REAL solutions:", real, "=> MK realisable over R:", bool(real))
# sanity: complex solution gives a genuine MK realisation with 8 distinct points (so the equations are right)
s0 = sols[0]
pts = {i: p[i].subs(s0) for i in range(8)}
distinct = all(sp.simplify(pts[i].cross(pts[j]).norm()) != 0 for i, j in itertools.combinations(range(8), 2))
inc = all(sp.simplify(sp.Matrix.hstack(*[pts[i] for i in sorted(L)]).det()) == 0 for L in MK)
print("complex solution: distinct points:", distinct, " all 8 incidences hold:", inc)
