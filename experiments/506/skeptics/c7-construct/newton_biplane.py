"""Try to REALISE the (7,4,2) biplane (Fano-line complements as concyclic/collinear quadruples)
in the real Mobius plane, numerically, with random restarts + Levenberg-Marquardt.
Normalisation: point 2 -> infinity, point 0 -> 0, point 1 -> 1 (Mobius group used up).
Blocks (must be concyclic or collinear): [0,1,2,3],[0,1,4,5],[2,3,4,5],[0,2,4,6],[1,3,4,6],[1,2,5,6],[0,3,5,6]
Blocks containing 2 (=inf) are straight lines through the other three points.
Unknowns: z3,z4,z5,z6 (8 reals). Residuals: 4 collinearity determinants + 3 concyclicity determinants,
normalised to be scale-free. Also a penalty to keep points apart is NOT used; instead we filter
solutions by pairwise distances afterwards."""
import numpy as np, sys
from scipy.optimize import least_squares
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
def coll(a,b,c):
    return ((b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1]))
def concyc(a,b,c,d):
    M=np.array([[p[0],p[1],p[0]**2+p[1]**2,1.0] for p in (a,b,c,d)])
    return np.linalg.det(M)
def pts(v):
    return {0:(0.0,0.0),1:(1.0,0.0),3:(v[0],v[1]),4:(v[2],v[3]),5:(v[4],v[5]),6:(v[6],v[7])}
def resid(v):
    P=pts(v)
    r=[coll(P[0],P[1],P[3]), coll(P[3],P[4],P[5]), coll(P[0],P[4],P[6]), coll(P[1],P[5],P[6]),
       concyc(P[0],P[1],P[4],P[5]), concyc(P[1],P[3],P[4],P[6]), concyc(P[0],P[3],P[5],P[6])]
    return np.array(r)
def mindist(v):
    P=pts(v); ks=list(P)
    return min(np.hypot(P[a][0]-P[b][0],P[a][1]-P[b][1]) for i,a in enumerate(ks) for b in ks[i+1:])
def maxnorm(v):
    P=pts(v); return max(np.hypot(*P[k]) for k in P)
best=[]
N=int(sys.argv[2]) if len(sys.argv)>2 else 3000
for it in range(N):
    v0=rng.normal(size=8)*rng.choice([0.5,1,3,10])
    try:
        sol=least_squares(resid,v0,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=4000)
    except Exception as e:
        continue
    v=sol.x; r=np.abs(resid(v)).max(); md=mindist(v); mx=maxnorm(v)
    # scale-free residual: divide by size^4 for concyclicity, size^2 for collinearity
    s=max(mx,1.0)
    rr=np.abs(resid(v))/np.array([s**2]*4+[s**4]*3)
    best.append((rr.max(), md/s, v))
best.sort(key=lambda t:t[0])
print("top 15 by scale-free residual (residual, mindist/size):")
for r,md,v in best[:15]:
    print(f"  res={r:.3e}  mindist/size={md:.3e}  v={np.round(v,4)}")
good=[(r,md,v) for r,md,v in best if r<1e-9 and md>1e-4]
print("non-degenerate solutions (res<1e-9, mindist/size>1e-4):",len(good))
for r,md,v in good[:10]: print("  ",r,md,np.round(v,6))
