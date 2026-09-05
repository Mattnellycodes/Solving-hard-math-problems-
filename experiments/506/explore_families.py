"""Exploration for Erdos #506: few circles determined by n points (not all on a circle/line).
Mobius reformulation: let B(P) = set of 'blocks' = circles or lines through >=3 points of P.
For an inversion centre O not in P, the Euclidean count of the inverted set is |B(P)| - deg(O),
where deg(O) = number of blocks through O (blocks through O become lines). With O = infinity we get
the Euclidean count of P itself, deg(inf) = #lines. So best(P) = |B(P)| - max_O deg(O), where O ranges
over P-free points; the maximum is attained at pairwise intersections of blocks (or infinity)."""
import math, itertools, sys, random
from collections import defaultdict
EPS=1e-7
def key(v): return round(v/EPS)
def block(a,b,c):
    ax,ay=a; bx,by=b; cx,cy=c
    d=2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    if abs(d)<1e-9:
        # line: normalize (A,B,C) with A x + B y = C
        A=by-ay; B=ax-bx; C=A*ax+B*ay
        nrm=math.hypot(A,B); A/=nrm; B/=nrm; C/=nrm
        if A<-1e-12 or (abs(A)<1e-12 and B<0): A,B,C=-A,-B,-C
        return ('L',key(A),key(B),key(C))
    a2=ax*ax+ay*ay; b2=bx*bx+by*by; c2=cx*cx+cy*cy
    ux=(a2*(by-cy)+b2*(cy-ay)+c2*(ay-by))/d
    uy=(a2*(cx-bx)+b2*(ax-cx)+c2*(bx-ax))/d
    r=math.hypot(ax-ux,ay-uy)
    return ('C',key(ux),key(uy),key(r))
def blocks(pts):
    Bl=defaultdict(set)
    for i,j,k in itertools.combinations(range(len(pts)),3):
        Bl[block(pts[i],pts[j],pts[k])].update((i,j,k))
    return Bl
def circ_line_intersections(b1,b2):
    """intersection points of two blocks (floats)."""
    out=[]
    if b1[0]=='L' and b2[0]=='L':
        A1,B1,C1=[x*EPS for x in b1[1:]]; A2,B2,C2=[x*EPS for x in b2[1:]]
        det=A1*B2-A2*B1
        if abs(det)<1e-9: return out
        out.append(((C1*B2-C2*B1)/det,(A1*C2-A2*C1)/det))
        return out
    if b1[0]=='L': b1,b2=b2,b1
    cx,cy,r=[x*EPS for x in b1[1:]]
    if b2[0]=='L':
        A,B,C=[x*EPS for x in b2[1:]]
        d=A*cx+B*cy-C  # signed distance (A,B unit)
        if abs(d)>r+1e-9: return out
        h2=r*r-d*d
        if h2<0: h2=0
        h=math.sqrt(h2)
        px,py=cx-A*d, cy-B*d
        out.append((px+(-B)*h,py+A*h)); out.append((px-(-B)*h,py-A*h))
        return out
    dx,dy,r2=[x*EPS for x in b2[1:]]
    D=math.hypot(dx-cx,dy-cy)
    if D<1e-9 or D>r+r2+1e-9 or D<abs(r-r2)-1e-9: return out
    a=(r*r-r2*r2+D*D)/(2*D); h2=r*r-a*a
    if h2<0: h2=0
    h=math.sqrt(h2)
    mx,my=cx+a*(dx-cx)/D, cy+a*(dy-cy)/D
    out.append((mx+h*(dy-cy)/D, my-h*(dx-cx)/D)); out.append((mx-h*(dy-cy)/D, my+h*(dx-cx)/D))
    return out
def evaluate(pts):
    """returns (best_count, euclid_count, nblocks, best_O) ; None if degenerate"""
    n=len(pts)
    for p,q in itertools.combinations(pts,2):
        if abs(p[0]-q[0])<1e-9 and abs(p[1]-q[1])<1e-9: return None
    Bl=blocks(pts)
    if any(len(s)==n for s in Bl.values()): return None
    nb=len(Bl)
    lines=sum(1 for b in Bl if b[0]=='L')
    euclid=nb-lines
    # candidate O: intersections of pairs of blocks; count blocks through each
    bl=list(Bl.keys())
    cand=defaultdict(int)
    ptkeys=set((key(x),key(y)) for x,y in pts)
    seen=defaultdict(set)
    for i in range(len(bl)):
        for j in range(i+1,len(bl)):
            for (x,y) in circ_line_intersections(bl[i],bl[j]):
                k=(round(x/1e-5),round(y/1e-5))
                if (key(x),key(y)) in ptkeys: continue
                seen[k].add(i); seen[k].add(j)
    bestdeg=lines; bestO='inf'
    for k,s in seen.items():
        # verify membership of point k in each block robustly
        x,y=k[0]*1e-5,k[1]*1e-5
        deg=0
        for b in bl:
            if b[0]=='L':
                A,B,C=[v*EPS for v in b[1:]]
                if abs(A*x+B*y-C)<1e-5: deg+=1
            else:
                cx,cy,r=[v*EPS for v in b[1:]]
                if abs(math.hypot(x-cx,y-cy)-r)<1e-5: deg+=1
        if deg>bestdeg: bestdeg=deg; bestO=(x,y)
    return (nb-bestdeg, euclid, nb, bestO)
def formula(n): return (n-1)*(n-2)//2+1-(n-1)//2
best={}
def consider(name,pts,verbose=False):
    n=len(pts)
    if n<4 or n>20: return
    r=evaluate(pts)
    if r is None: return
    c,eu,nb,O=r
    if n not in best or c<best[n][0]:
        best[n]=(c,eu,nb,O,name)
        if verbose: print("new best",n,c,name,flush=True)
def poly(m,r,rot,cx=0,cy=0):
    return [(cx+r*math.cos(2*math.pi*i/m+rot), cy+r*math.sin(2*math.pi*i/m+rot)) for i in range(m)]
def seg_inter(p1,p2,p3,p4):
    x1,y1=p1;x2,y2=p2;x3,y3=p3;x4,y4=p4
    d=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
    if abs(d)<1e-12: return None
    t=((x1-x3)*(y3-y4)-(y1-y3)*(x3-x4))/d
    return (x1+t*(x2-x1), y1+t*(y2-y1))
def dedupe(pts):
    out=[]
    for p in pts:
        if all(abs(p[0]-q[0])>1e-7 or abs(p[1]-q[1])>1e-7 for q in out): out.append(p)
    return out
ratios=[2,3,1.5,math.sqrt(2),math.sqrt(3),(1+math.sqrt(5))/2,2.5,1/math.cos(math.pi/4),1/math.cos(math.pi/3),1/math.cos(math.pi/5),1/math.cos(math.pi/6),1/math.cos(math.pi/8),1/math.cos(math.pi/10),1/math.cos(math.pi/12)]
# (A) concentric polygons
for m in range(3,11):
    for k in (1,2,3):
        for rots in itertools.product([0,1],repeat=k):
            for rs in itertools.product(range(len(ratios)),repeat=k-1):
                radii=[1.0]
                for i in rs: radii.append(radii[-1]*ratios[i])
                pts=[]
                for layer in range(k): pts+=poly(m,radii[layer],rots[layer]*math.pi/m)
                for centre in (0,1):
                    P=pts+([(0.0,0.0)] if centre else [])
                    consider(f"{k}x{m}-gon rot={rots} radii={[round(x,3) for x in radii]} centre={centre}",P)
print("done A",flush=True)
# (B) star polygons: m-gon + intersections of step-s diagonals
for m in range(4,11):
    V=poly(m,1.0,0)
    for s in range(2,m//2+1):
        diag=[(V[i],V[(i+s)%m]) for i in range(m)]
        inter=[]
        for d1,d2 in itertools.combinations(diag,2):
            q=seg_inter(*d1,*d2)
            if q is not None and math.hypot(*q)<0.999: inter.append(q)
        inter=dedupe(inter)
        for centre in (0,1):
            P=dedupe(V+inter+([(0.0,0.0)] if centre else []))
            consider(f"star {m}/{s} +inner({len(inter)}) centre={centre}",P)
            P2=dedupe(inter+([(0.0,0.0)] if centre else []))
            consider(f"star {m}/{s} inner only({len(inter)}) centre={centre}",P2)
print("done B",flush=True)
# (C) grids and cube/hypercube projections
for a in range(2,6):
    for b in range(2,6):
        P=[(float(i),float(j)) for i in range(a) for j in range(b)]
        consider(f"grid {a}x{b}",P)
        consider(f"grid {a}x{b} + centre",dedupe(P+[((a-1)/2,(b-1)/2)]))
for d in [(1,1,1),(1,1,0.5),(1,0.5,0.3),(1,2,3),(0.2,0.3,1),(1,1,0.001),(1,0.3,0)]:
    nd=math.sqrt(sum(y*y for y in d)); d=[x/nd for x in d]
    u=[-d[1],d[0],0]; nu=math.hypot(*u)
    if nu<1e-9: continue
    u=[x/nu for x in u]
    v=[d[1]*u[2]-d[2]*u[1], d[2]*u[0]-d[0]*u[2], d[0]*u[1]-d[1]*u[0]]
    for dim in (3,4):
        P=[]
        for c in itertools.product([0,1],repeat=dim):
            cc=list(c)+[0]*(3-dim) if dim<3 else c
            # for dim 4 use a 4d->3d projection first: (x,y,z,w)->(x+0.3w, y+0.6w, z+0.9w)
            if dim==4: cc=(c[0]+0.3*c[3], c[1]+0.6*c[3], c[2]+0.9*c[3])
            P.append((sum(cc[i]*u[i] for i in range(3)), sum(cc[i]*v[i] for i in range(3))))
        consider(f"{dim}-cube proj d={[round(x,3) for x in d]}",dedupe(P))
print("done C",flush=True)
# (D) greedy deletion from rich seeds
seeds=[("grid4x4",[(float(i),float(j)) for i in range(4) for j in range(4)]),
       ("grid5x5",[(float(i),float(j)) for i in range(5) for j in range(5)]),
       ("2x8gon+c",poly(8,1,0)+poly(8,math.sqrt(2),0)+[(0.0,0.0)]),
       ("2x6gon+c",poly(6,1,0)+poly(6,2,0)+[(0.0,0.0)]),
       ("3x4gon+c",poly(4,1,0)+poly(4,2,0)+poly(4,3,0)+[(0.0,0.0)]),
       ("2x5+2x5",poly(5,1,0)+poly(5,(1+math.sqrt(5))/2,math.pi/5)+poly(5,2.2,0)+poly(5,2.2*(1+math.sqrt(5))/2,math.pi/5)),
       ("2x10gon",poly(10,1,0)+poly(10,1/math.cos(math.pi/10),math.pi/10)),
       ("grid4x4+diagcentres",[(float(i),float(j)) for i in range(4) for j in range(4)]+[(1.5,1.5)]),
       ]
for name,S in seeds:
    S=dedupe(S)
    cur=list(S)
    consider(name,cur)
    while len(cur)>5:
        bestc=None
        for i in range(len(cur)):
            Q=cur[:i]+cur[i+1:]
            r=evaluate(Q)
            if r is None: continue
            if bestc is None or r[0]<bestc[0]: bestc=(r[0],i)
        if bestc is None: break
        cur=cur[:bestc[1]]+cur[bestc[1]+1:]
        consider(f"{name} greedy-del -> {len(cur)}",cur)
print("done D",flush=True)
for n in sorted(best):
    c,eu,nb,O,name=best[n]
    flag="  <-- BELOW FORMULA" if c<formula(n) else ""
    print(f"n={n:2d} formula={formula(n):3d} best={c:3d} (euclid={eu}, blocks={nb}, O={O if O=='inf' else tuple(round(v,4) for v in O)}) {name}{flag}")
