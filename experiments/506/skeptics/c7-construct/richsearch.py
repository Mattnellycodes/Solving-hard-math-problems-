"""Search 7-subsets of many 'rich' planar configurations (float, tolerance-based counting) for
few circles. Anything with <= 11 circles is printed; <= 10 would refute the claim."""
import numpy as np, itertools, math, sys
TOL=1e-7
def circ(a,b,c):
    (x1,y1),(x2,y2),(x3,y3)=a,b,c
    d=2*((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1))
    sc=max(1.0,abs(x1),abs(y1),abs(x2),abs(y2),abs(x3),abs(y3))
    if abs(d)<1e-9*sc*sc: return None
    s1=x1*x1+y1*y1; s2=x2*x2+y2*y2; s3=x3*x3+y3*y3
    ux=((s2-s1)*(y3-y1)-(s3-s1)*(y2-y1))/d; uy=((x2-x1)*(s3-s1)-(x3-x1)*(s2-s1))/d
    r=math.hypot(x1-ux,y1-uy)
    return (ux,uy,r)
def count(P):
    keys=set(); nl=0
    for a,b,c in itertools.combinations(P,3):
        k=circ(a,b,c)
        if k is None: nl+=1; continue
        keys.add((round(k[0]/TOL),round(k[1]/TOL),round(k[2]/TOL)))
    return len(keys),nl
def degenerate(P):
    n,_=count(P)
    if n==0: return True
    if n>1: return False
    # single circle: all points on it?
    for a,b,c in itertools.combinations(P,3):
        k=circ(a,b,c)
        if k: return all(abs(math.hypot(x-k[0],y-k[1])-k[2])<1e-6 for x,y in P)
    return True
def invert(P,c,R2=1.0):
    out=[]
    for x,y in P:
        dx,dy=x-c[0],y-c[1]; d2=dx*dx+dy*dy
        if d2<1e-12: continue
        out.append((c[0]+R2*dx/d2,c[1]+R2*dy/d2))
    return out
def dedup(P):
    out=[]
    for p in P:
        if all(math.hypot(p[0]-q[0],p[1]-q[1])>1e-9 for q in out): out.append(p)
    return out
def search(name,P,k=7,report=11):
    P=dedup(P); best=(999,None)
    if len(P)<k: return best
    for S in itertools.combinations(P,k):
        n,nl=count(S)
        if n<=report and not degenerate(S):
            if n<best[0]: best=(n,S)
            if n<=10: print("!!! <=10:",name,n,S,flush=True)
    print(f"{name:45s} N={len(P):2d} best 7-subset circles={best[0]}", "" if best[1] is None else np.round(best[1],4).tolist(), flush=True)
    return best
fams={}
def polygon(m,r=1.0,rot=0.0,c=(0,0)):
    return [(c[0]+r*math.cos(rot+2*math.pi*i/m),c[1]+r*math.sin(rot+2*math.pi*i/m)) for i in range(m)]
for m in range(5,13):
    fams[f"{m}-gon+center"]=polygon(m)+[(0,0)]
    fams[f"{m}-gon+center+midpoints"]=polygon(m)+[(0,0)]+polygon(m,math.cos(math.pi/m),math.pi/m)
    fams[f"{m}-gon+concentric {m}-gon r=2"]=polygon(m)+polygon(m,2.0)+[(0,0)]
    fams[f"{m}-gon+concentric rotated"]=polygon(m)+polygon(m,2.0,math.pi/m)+[(0,0)]
# triangle centres family
def triangle_family(A,B,C):
    A,B,C=map(np.array,(A,B,C))
    def foot(P,Q,R):  # foot of perpendicular from P to QR
        d=R-Q; t=np.dot(P-Q,d)/np.dot(d,d); return Q+t*d
    D=foot(A,B,C);E=foot(B,C,A);F=foot(C,A,B)
    # orthocentre
    H=A+B+C-2*np.array([0,0]) # placeholder
    # circumcentre
    ax,ay=A;bx,by=B;cx,cy=C
    d=2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    ux=((ax*ax+ay*ay)*(by-cy)+(bx*bx+by*by)*(cy-ay)+(cx*cx+cy*cy)*(ay-by))/d
    uy=((ax*ax+ay*ay)*(cx-bx)+(bx*bx+by*by)*(ax-cx)+(cx*cx+cy*cy)*(bx-ax))/d
    O=np.array([ux,uy]); H=A+B+C-2*O; N=(O+H)/2; G=(A+B+C)/3
    Ma,Mb,Mc=(B+C)/2,(C+A)/2,(A+B)/2; Ea,Eb,Ec=(A+H)/2,(B+H)/2,(C+H)/2
    return [tuple(p) for p in (A,B,C,H,O,N,G,D,E,F,Ma,Mb,Mc,Ea,Eb,Ec)]
fams["triangle(0,3),(1,0),(3,0) centres"]=triangle_family((0,3),(1,0),(3,0))
fams["triangle(0,0),(20,0),(5,15) centres"]=triangle_family((0,0),(20,0),(5,15))
fams["triangle(0,0),(4,0),(1,3) centres"]=triangle_family((0,0),(4,0),(1,3))
fams["triangle(0,0),(6,0),(2,5) centres"]=triangle_family((0,0),(6,0),(2,5))
fams["3x3 grid"]=[(x,y) for x in range(3) for y in range(3)]
fams["4x4 grid"]=[(x,y) for x in range(4) for y in range(4)]
fams["two squares+centre+4"]=[(1,0),(0,1),(-1,0),(0,-1),(2,0),(0,2),(-2,0),(0,-2),(0,0),(1,1),(-1,1),(1,-1),(-1,-1)]
fams["record8 + extras"]=[(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10),(0,7.5),(3,7.5),(5,7.5),(1.5,7.5)]
fams["n7 construction + extras"]=[(0,0),(20,0),(5,15),(5,0),(2,6),(10,10),(5,5),(10,0),(12.5,7.5),(2.5,7.5),(5,10),(7.5,2.5),(3.5,10.5)]
# hexagonal lattice patch
fams["hex lattice 12"]=[(i+0.5*j, j*math.sqrt(3)/2) for i in range(-2,3) for j in range(-1,2) if abs(i+0.5*j)<2.6][:14]
# cube / octahedron / cuboctahedron projections (stereographic from a generic centre) -- Segre style
import random
random.seed(1)
def stereo(V,c=None):
    V=[np.array(v,float)/np.linalg.norm(v) for v in V]
    if c is None: c=np.array([0.3,0.2,0.93]); c/=np.linalg.norm(c)
    # rotate so c -> north pole then project from north pole
    z=c; x=np.cross(z,[0,0,1.0]); 
    if np.linalg.norm(x)<1e-9: x=np.array([1.0,0,0])
    x/=np.linalg.norm(x); y=np.cross(z,x)
    out=[]
    for v in V:
        vx,vy,vz=np.dot(v,x),np.dot(v,y),np.dot(v,z)
        if abs(1-vz)<1e-9: continue
        out.append((vx/(1-vz),vy/(1-vz)))
    return out
cube=[(sx,sy,sz) for sx in (-1,1) for sy in (-1,1) for sz in (-1,1)]
octa=[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
cubocta=[(a,b,0) for a in (-1,1) for b in (-1,1)]+[(a,0,b) for a in (-1,1) for b in (-1,1)]+[(0,a,b) for a in (-1,1) for b in (-1,1)]
phi=(1+5**.5)/2
ico=[(0,s1,s2*phi) for s1 in (-1,1) for s2 in (-1,1)]+[(s1,s2*phi,0) for s1 in (-1,1) for s2 in (-1,1)]+[(s2*phi,0,s1) for s1 in (-1,1) for s2 in (-1,1)]
for nm,V in [("cube",cube),("octa",octa),("cubocta",cubocta),("ico",ico),("cube+octa",cube+[tuple(3**.5*np.array(o)) for o in octa])]:
    for k in range(3):
        c=np.random.default_rng(k).normal(size=3)
        fams[f"stereo {nm} c{k}"]=stereo(V,c)
    # project from a vertex direction (vertex -> infinity)
    fams[f"stereo {nm} from vertex"]=stereo(V,np.array(V[0],float))
results={}
for name,P in fams.items():
    P=[(float(x),float(y)) for x,y in P]
    results[name]=search(name,P)
    # inversions about the centroid and about each of first 3 points
    cen=(sum(x for x,y in P)/len(P),sum(y for x,y in P)/len(P))
    centers=[cen]+P[:3]+[(cen[0]+0.37,cen[1]-0.21)]
    for i,c in enumerate(centers):
        Q=invert(P,c,1.0+0.1*i)
        search(name+f" inv{i}",Q)
