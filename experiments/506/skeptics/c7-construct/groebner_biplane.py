"""Exact algebraic check (own parametrisation, independent of the theory agent's).
Complete quadrilateral: L_a: y=0; L_c: y=mc*x (through 0); L_d: y=md*(x-1) (through 1);
L_b: y=mb*(x-t) (through (t,0)).  Vertices: P0=(0,0)=La∩Lc, P1=(1,0)=La∩Ld, P3=(t,0)=La∩Lb,
P4=Lb∩Lc, P5=Lb∩Ld, P6=Lc∩Ld.  Diagonal pairs (0,5),(1,4),(3,6).
Require {0,1,4,5},{1,3,4,6},{0,3,5,6} concyclic (determinant = 0, incl. collinear).
This is exactly the structure at a point of the biplane (point 2 sent to infinity).
Rotation freedom WAS used to make L_a horizontal, translation+scaling to put P0=(0,0),P1=(1,0);
so no generality is lost except: lines parallel to L_a would have slope 0 (excluded since each
pair of lines meets at a vertex), vertical lines cannot be written y=m x + c -> handled separately."""
import sympy as sp
mb,mc,md,t,u=sp.symbols('mb mc md t u')
def inter(m1,c1,m2,c2):  # y=m1 x + c1 , y = m2 x + c2
    x=(c2-c1)/(m1-m2); return (sp.simplify(x), sp.simplify(m1*x+c1))
P={0:(sp.Integer(0),sp.Integer(0)),1:(sp.Integer(1),sp.Integer(0)),3:(t,sp.Integer(0))}
P[4]=inter(mb,-mb*t,mc,0); P[5]=inter(mb,-mb*t,md,-md); P[6]=inter(mc,0,md,-md)
def cyc(ids):
    M=sp.Matrix([[P[i][0],P[i][1],P[i][0]**2+P[i][1]**2,1] for i in ids])
    return sp.factor(sp.numer(sp.together(M.det())))
Q=[cyc([0,1,4,5]),cyc([1,3,4,6]),cyc([0,3,5,6])]
for q in Q: print("factored:",q)
# core factors: strip the obvious degeneracy factors
def core(q):
    fs=sp.factor_list(q)[1]
    keep=[]
    for f,e in fs:
        if f in (mb,mc,md,t,t-1,mb-mc,mb-md,mc-md) or (-f) in (mb,mc,md,t,t-1,mb-mc,mb-md,mc-md): continue
        keep.append(f**e)
    return sp.Mul(*keep)
C=[sp.expand(core(q)) for q in Q]
print("core polynomials:")
for c in C: print("  ",sp.factor(c))
G=sp.groebner(C,mb,mc,md,t,order='lex')
print("Groebner (lex) of the three core polynomials, NO saturation:")
for g in G: print("  ",sp.factor(g))
# solutions
sols=sp.solve(C,[mb,mc,md,t],dict=True)
print("sympy solve components:",len(sols))
for s in sols: print("  ",s)
# what do the components mean? check degeneracy conditions on each
print("vertex coords:",{k:P[k] for k in (4,5,6)})
