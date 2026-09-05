import sys, math, itertools
sys.path.insert(0,'/home/user/Solving-hard-math-problems-/experiments/506')
import explore_families as ef, circles_exact as ce
s3=math.sqrt(3)
# image of cube under stereographic projection from (0,0,sqrt3): squares of half-side (3-+sqrt3)/2
a=(3+s3)/2; b=(3-s3)/2
img=[(sx*a,sy*a) for sx in (-1,1) for sy in (-1,1)]+[(sx*b,sy*b) for sx in (-1,1) for sy in (-1,1)]
print("cube image (pole on 4-fold axis), evaluate:",ef.evaluate(img))
# same with rotated-by-45deg ratio 2 configuration (our record)
ours=[(1,0),(-1,0),(0,1),(0,-1),(2,0),(-2,0),(0,2),(0,-2)]
print("two squares ratio 2 (ours), evaluate:",ef.evaluate([(float(x),float(y)) for x,y in ours]))
print("two squares ratio 2 exact count:",ce.count_circles(ours), "degenerate?",ce.all_on_one_circle_or_line(ours))
# ratio r family: squares (+-1,+-1) and (+-r,+-r) -- exact rational check for several r
for r in [2,3,5,7,'5/2','7/3','10/3']:
    from fractions import Fraction as Fr
    R=Fr(r)
    pts=[(sx,sy) for sx in (-1,1) for sy in (-1,1)]+[(sx*R,sy*R) for sx in (-1,1) for sy in (-1,1)]
    print(" ratio",r,"exact circles:",ce.count_circles(pts))
# generic pole projection of cube
def stereo(pole,v):
    N=pole; t=1/(1-sum(N[i]*v[i] for i in range(3))/3)
    w=[N[i]+t*(v[i]-N[i]) for i in range(3)]
    u=[N[i]/s3 for i in range(3)]
    e1=[-u[1],u[0],0]; n=math.hypot(*e1); e1=[c/n for c in e1]
    e2=[u[1]*e1[2]-u[2]*e1[1], u[2]*e1[0]-u[0]*e1[2], u[0]*e1[1]-u[1]*e1[0]]
    return (sum(w[i]*e1[i] for i in range(3)), sum(w[i]*e2[i] for i in range(3)))
V=list(itertools.product([-1,1],repeat=3))
th,ph=0.7,1.3
pole=[s3*math.sin(th)*math.cos(ph), s3*math.sin(th)*math.sin(ph), s3*math.cos(th)]
print("generic pole:",ef.evaluate([stereo(pole,v) for v in V]))
# pole on a face circle only (z=1 plane, non-vertex): (sqrt2, 0, 1)
pole=[math.sqrt(2),0.0,1.0]
print("pole on one face-circle:",ef.evaluate([stereo(pole,v) for v in V]))
