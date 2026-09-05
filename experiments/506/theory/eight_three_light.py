"""Light, fully explicit proof that no 8 real points have 8 three-point lines (t_3(8) <= 7).
Step 1 (combinatorial, brute force): every family of 8 triples on 8 points pairwise sharing <= 1
point is isomorphic to the Mobius-Kantor configuration MK = {i, i+1, i+3 mod 8}.
Step 2 (algebraic, exact): a projective realisation of MK with the 4-point frame {0,1,2,5} (no three
collinear in MK) forces  p3=(1,1,0), p6=(a,1,1), p7=(b,0,1), p4=(0,c,1) with a=b=1-c and c^2-c+1=0,
which has no real root.  Frame points and chart choices: p6 lies on line(p5,p0)={y=z}: its points are
(a,1,1) or p0 itself; p7 on line(p0,p2)={y=0}: (b,0,1) or p0; p4 on line(p1,p2)={x=0}: (0,c,1) or p1;
p3 = line(p0,p1) ∩ line(p2,p5) = (1,1,0).  Points must be distinct, so the alternatives are excluded."""
import itertools, sympy as sp
n=8
triples=[frozenset(t) for t in itertools.combinations(range(n),3)]
perms=list(itertools.permutations(range(n)))
def canon(F): return min(tuple(sorted(tuple(sorted(s[i] for i in t)) for t in F)) for s in perms)
classes=set()
def rec(F,start):
    if len(F)==8: classes.add(canon(F)); return
    for j in range(start,len(triples)):
        t=triples[j]
        if all(len(t&u)<=1 for u in F):
            if all(sum(1 for u in F if p in u)<3 for p in t): rec(F+[t],j+1)
rec([triples[0]],1)
MK=[frozenset({i,(i+1)%8,(i+3)%8}) for i in range(8)]
print("iso classes of 8-triple systems on 8 points (pairwise <=1):",len(classes),
      "; Mobius-Kantor is the class:", canon(MK) in classes)
a,b,c=sp.symbols('a b c')
p={0:sp.Matrix([1,0,0]),1:sp.Matrix([0,1,0]),2:sp.Matrix([0,0,1]),5:sp.Matrix([1,1,1]),
   3:sp.Matrix([1,1,0]),6:sp.Matrix([a,1,1]),7:sp.Matrix([b,0,1]),4:sp.Matrix([0,c,1])}
eqs=[sp.Matrix.hstack(*[p[i] for i in sorted(L)]).det() for L in MK]
eqs=[sp.expand(e) for e in eqs]
print("incidence equations:",[e for e in eqs if e!=0])
G=sp.groebner([e for e in eqs if e!=0],a,b,c,order='lex')
print("Groebner basis (lex a>b>c):",list(G))
poly=[g for g in G if g.free_symbols=={c}][0]
print("univariate polynomial in c:",sp.factor(poly),"  real roots:",sp.real_roots(sp.Poly(poly,c)))
