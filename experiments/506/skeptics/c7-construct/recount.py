"""Independent exact recount (Fractions) of circles through >=3 points. Written from scratch.
A circle is identified by the exact (center, r^2). Collinear triples are skipped."""
from fractions import Fraction as F
from itertools import combinations
import sys

def circ(a,b,c):
    (x1,y1),(x2,y2),(x3,y3)=a,b,c
    d = 2*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
    if d == 0: return None
    s1=x1*x1+y1*y1; s2=x2*x2+y2*y2; s3=x3*x3+y3*y3
    ux = ((s2-s1)*(y3-y1)-(s3-s1)*(y2-y1))/d
    uy = ((x2-x1)*(s3-s1)-(x3-x1)*(s2-s1))/d
    return (ux,uy,(x1-ux)**2+(y1-uy)**2)

def lines_through(pts):
    L={}
    for i,j,k in combinations(range(len(pts)),3):
        (x1,y1),(x2,y2),(x3,y3)=pts[i],pts[j],pts[k]
        if (x2-x1)*(y3-y1)-(x3-x1)*(y2-y1)==0:
            # canonical line: (a,b,c) with a x + b y = c, normalised
            a=y2-y1; b=x1-x2; c=a*x1+b*y1
            if a!=0: b,c,a=b/a,c/a,F(1)
            else: c,b=c/b,F(1)
            L.setdefault((a,b,c),set()).update((i,j,k))
    return L

def report(pts):
    pts=[(F(x),F(y)) for x,y in pts]
    assert len(set(pts))==len(pts)
    C={}
    for i,j,k in combinations(range(len(pts)),3):
        c=circ(pts[i],pts[j],pts[k])
        if c is not None: C.setdefault(c,set()).update((i,j,k))
    L=lines_through(pts)
    degenerate = (len(C)==0) or any(len(s)==len(pts) for s in C.values())
    return len(C), sorted(len(s) for s in C.values()), sorted(len(s) for s in L.values()), degenerate

if __name__=="__main__":
    P6=[(0,0),(20,0),(5,15),(5,0),(2,6),(10,10)]
    P7=P6+[(5,5)]
    P8=[(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)]
    for name,P in [("n6",P6),("n7",P7),("n8",P8)]:
        n,cs,ls,deg=report(P)
        print(name, "circles =",n, "circle sizes",cs, "line sizes",ls, "degenerate?",deg)
