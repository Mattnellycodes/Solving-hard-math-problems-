"""Own proof that the Möbius-Kantor configuration (8_3) = {i, i+1, i+3 mod 8} has no realisation by
8 distinct points of RP^2 whose only collinear triples are the 8 configuration lines (the only kind of
realisation that can occur as the 3-point-line structure of a real point set).
Frame: p0=(1,0,0), p1=(0,1,0), p2=(0,0,1), p5=(1,1,1): no three of {0,1,2,5} lie on a configuration
line, hence they are in general position.  Lines: {0,1,3}: z=0; {2,3,5}: x=y  => p3=(1,1,0).
{5,6,0}: y=z => p6=(a,1,1) (y=0 would give p6=p0); {7,0,2}: y=0 => p7=(b,0,1) (z=0 gives p7=p0);
{1,2,4}: x=0 => p4=(0,c,1) (z=0 gives p4=p1).  Remaining lines {3,4,6}, {4,5,7}, {6,7,1}."""
import sympy as sp
from itertools import combinations
a, b, c = sp.symbols("a b c")
P = {0: (1, 0, 0), 1: (0, 1, 0), 2: (0, 0, 1), 5: (1, 1, 1), 3: (1, 1, 0), 6: (a, 1, 1), 7: (b, 0, 1), 4: (0, c, 1)}
lines = [((i) % 8, (i + 1) % 8, (i + 3) % 8) for i in range(8)]
det = lambda tr: sp.expand(sp.Matrix([P[p] for p in tr]).det())
used = [(0, 1, 3), (2, 3, 5), (5, 6, 0), (7, 0, 2), (1, 2, 4)]
for tr in used:
    assert det(tr) == 0, tr
rest = [tr for tr in lines if tr not in used]
print("remaining line conditions:", {tr: det(tr) for tr in rest})
G = sp.groebner([det(tr) for tr in rest], a, b, c, order="lex")
print("lex Gröbner basis:", G.exprs)
last = G.exprs[-1]
print("univariate polynomial in c:", last, " real roots:", sp.real_roots(sp.Poly(last, c)))
assert not sp.real_roots(sp.Poly(last, c))
print("(8_3) has no real realisation of the required kind")
