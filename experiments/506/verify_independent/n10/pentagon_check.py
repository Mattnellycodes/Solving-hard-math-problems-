"""Main agent's numerical sanity check (2026-09-05) of the decisive geometric step for c(10) >= 33.
The only abstract structure surviving all combinatorial filters is: two disjoint 5-blocks S={0..4}, R={5..9},
the 20 four-blocks below, and 10 three-point lines. The agents claim (a) every realisation of the 22 rich
blocks is two concentric regular pentagons (aligned or rotated by pi/5) with a free radius ratio rho, and
(b) for no rho does a point O exist lying on the 10 circumcircles of the line-triples (O would be the
inversion centre sending those circles to lines). We check (a) by computing the block structure of concentric
pentagons for sample rho and matching it to the abstract family, and (b) by scanning rho finely: intersect
the first two circumcircles and measure the distance of the intersection points to the other eight circles.
"""
import math, itertools, numpy as np
from collections import defaultdict
F4 = [[0,1,6,9],[0,1,7,8],[0,2,5,9],[0,2,6,8],[0,3,5,7],[0,3,8,9],[0,4,5,6],[0,4,7,9],[1,2,5,8],[1,2,7,9],[1,3,5,9],[1,3,6,7],[1,4,5,7],[1,4,6,8],[2,3,5,6],[2,3,7,8],[2,4,6,7],[2,4,8,9],[3,4,5,8],[3,4,6,9]]
LINES = [[0,1,5],[0,4,8],[0,6,7],[1,3,7],[1,8,9],[2,3,9],[2,4,6],[2,7,8],[3,5,6],[4,5,9]]
Fabs = [frozenset(b) for b in F4]
def circ(a,b,c):
    ax,ay=a; bx,by=b; cx,cy=c
    d=2*(ax*(by-cy)+bx*(cy-ay)+cx*(ay-by))
    if abs(d)<1e-12: return None
    a2=ax*ax+ay*ay; b2=bx*bx+by*by; c2=cx*cx+cy*cy
    ux=(a2*(by-cy)+b2*(cy-ay)+c2*(ay-by))/d; uy=(a2*(cx-bx)+b2*(ax-cx)+c2*(bx-ax))/d
    return (ux,uy,math.hypot(ax-ux,ay-uy))
def blocks(pts,tol=1e-7):
    B=defaultdict(set)
    for i,j,k in itertools.combinations(range(len(pts)),3):
        c=circ(pts[i],pts[j],pts[k])
        if c is None: key=('L',)
        else: key=('C',round(c[0]/tol),round(c[1]/tol),round(c[2]/tol))
        B[key].update((i,j,k))
    return [frozenset(s) for s in B.values() if len(s)>=4]
def pent(rho,phi):
    P=[(math.cos(2*math.pi*k/5),math.sin(2*math.pi*k/5)) for k in range(5)]
    Q=[(rho*math.cos(2*math.pi*k/5+phi),rho*math.sin(2*math.pi*k/5+phi)) for k in range(5)]
    return P+Q
def find_iso(geo4):
    """find labelings (outer->S, inner->R) mapping geometric 4-blocks onto Fabs"""
    G=set(geo4); isos=[]
    for p in itertools.permutations(range(5)):
        for q in itertools.permutations(range(5)):
            lab={k:p[k] for k in range(5)}; lab.update({5+k:5+q[k] for k in range(5)})
            img=set(frozenset(lab[x] for x in b) for b in G)
            if img==set(Fabs): isos.append(lab)
    return isos
for phi in (0.0, math.pi/5):
    for rho in (0.37, 0.61, 0.8):
        pts=pent(rho,phi); B=blocks(pts)
        four=[b for b in B if len(b)==4]; big=[b for b in B if len(b)>=5]
        print(f"phi={phi:.3f} rho={rho}: blocks>=5: {[len(b) for b in big]}, #4-blocks={len(four)}", end=' ')
        isos=find_iso(four) if len(four)==20 else []
        print(f"isomorphisms with abstract F: {len(isos)}")
        if isos:
            lab=isos[0]; inv={v:k for k,v in lab.items()}
            # (b) concurrency scan over rho for this labeling family
            best=(1e9,None)
            for r in np.linspace(0.02,0.999,4000):
                pp=pent(r,phi)
                tri=[[pp[inv[x]] for x in L] for L in LINES]
                cs=[circ(*t) for t in tri]
                if any(c is None for c in cs): continue
                # intersections of circles 0 and 1
                (x0,y0,r0),(x1,y1,r1)=cs[0],cs[1]; D=math.hypot(x1-x0,y1-y0)
                if D<1e-12 or D>r0+r1 or D<abs(r0-r1): continue
                a=(r0*r0-r1*r1+D*D)/(2*D); h=math.sqrt(max(r0*r0-a*a,0)); mx=x0+a*(x1-x0)/D; my=y0+a*(y1-y0)/D
                for sgn in (1,-1):
                    O=(mx+sgn*h*(y1-y0)/D, my-sgn*h*(x1-x0)/D)
                    err=max(abs(math.hypot(O[0]-c[0],O[1]-c[1])-c[2]) for c in cs[2:])
                    if err<best[0]: best=(err,r)
            print(f"   min over rho of max-distance of an intersection point of circles 1,2 to circles 3..10: {best[0]:.4g} at rho={best[1]:.4f}")
            # infinity as inversion centre: all 10 triples collinear?
            print("   (O = infinity would need all 10 line-triples collinear; they are not for any rho since points of a regular pentagon are never collinear with the other pentagon's points except by coincidence — checked at sample rho:)",
                  all(circ(*[pent(0.61,phi)[inv[x]] for x in L]) is not None for L in LINES))
