"""Exact check of Segre's cube example for Erdos #506.
Cube vertices (+-1,+-1,+-1) lie on the sphere x^2+y^2+z^2=3.  Circles on the sphere = planes
through >=3 vertices.  Stereographic projection from a pole N on the sphere (N not a vertex)
maps circles through N to lines, all other circles to circles.  So the Euclidean circle count of
the projected set is  (#planes through >=3 vertices) - (#such planes through N)."""
import itertools, sympy as sp
V=[sp.Matrix(v) for v in itertools.product([-1,1],repeat=3)]
# planes through >=3 vertices
planes={}
for a,b,c in itertools.combinations(range(8),3):
    n=(V[b]-V[a]).cross(V[c]-V[a])
    g=sp.gcd(list(n)); n=n/g
    # normalise sign
    for t in n:
        if t!=0:
            if t<0: n=-n
            break
    d=n.dot(V[a])
    k=(tuple(n),d)
    planes.setdefault(k,set()).update((a,b,c))
print("number of planes through >=3 vertices:",len(planes))
from collections import Counter
print("sizes:",Counter(len(s) for s in planes.values()))
# degree of non-vertex sphere points: intersect pairs of planes with the sphere
x,y,z=sp.symbols('x y z',real=True)
best=0; bestpts=[]
keys=list(planes)
def on_plane(p,k):
    n,d=k; return sp.simplify(sp.Matrix(n).dot(p)-d)==0
deg={}
for i,j in itertools.combinations(range(len(keys)),2):
    (n1,d1),(n2,d2)=keys[i],keys[j]
    sols=sp.solve([sp.Matrix(n1).dot(sp.Matrix([x,y,z]))-d1, sp.Matrix(n2).dot(sp.Matrix([x,y,z]))-d2, x**2+y**2+z**2-3],[x,y,z],dict=True)
    for s in sols:
        p=sp.Matrix([sp.nsimplify(s[x]),sp.nsimplify(s[y]),sp.nsimplify(s[z])])
        if any(p==v for v in V): continue
        key=tuple(sp.simplify(t) for t in p)
        if key in deg: continue
        deg[key]=sum(1 for k in keys if on_plane(p,k))
print("max degree of a non-vertex sphere point:",max(deg.values()))
print("points of max degree:",[k for k,v in deg.items() if v==max(deg.values())])
print("degree distribution:",Counter(deg.values()))
# stereographic projection from N=(0,0,sqrt3) onto z=0
N=sp.Matrix([0,0,sp.sqrt(3)])
def proj(v):
    t=N[2]/(N[2]-v[2])
    return (sp.simplify(N[0]+t*(v[0]-N[0])), sp.simplify(N[1]+t*(v[1]-N[1])))
img=[proj(v) for v in V]
print("image of cube under stereographic projection from (0,0,sqrt3):")
for v,p in zip(V,img): print("  ",tuple(v),"->",p)
import sys; sys.path.insert(0,'/home/user/Solving-hard-math-problems-/experiments/506')
import circles_exact as ce
try:
    print("exact analyze:",ce.analyze(img))
except Exception as e:
    print("circles_exact.analyze failed on sympy input:",e)
# float check with explore_families.evaluate
import explore_families as ef
pts=[(float(a),float(b)) for a,b in img]
print("evaluate (best_count, euclid_count, nblocks, best_O):",ef.evaluate(pts))
# generic pole
import random
random.seed(1)
th,ph=0.7,1.3
Ng=sp.Matrix([sp.sqrt(3)*sp.sin(th)*sp.cos(ph),sp.sqrt(3)*sp.sin(th)*sp.sin(ph),sp.sqrt(3)*sp.cos(th)])
def projg(v):
    # project from Ng onto plane through origin orthogonal to Ng, using orthonormal basis
    u=Ng/sp.sqrt(3); e1=sp.Matrix([-u[1],u[0],0]); e1=e1/e1.norm(); e2=u.cross(e1)
    t=1/(1-u.dot(v)/sp.sqrt(3))
    w=Ng+t*(v-Ng)
    return (float(w.dot(e1)),float(w.dot(e2)))
ptsg=[projg(v) for v in V]
print("generic pole evaluate:",ef.evaluate(ptsg))
