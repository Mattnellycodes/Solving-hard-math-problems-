# quick float Moebius evaluator in Python to compare with C
import itertools, math
def blocks(P):
    n=len(P); bl=[]
    done=[]
    for i,j,k in itertools.combinations(range(n),3):
        mk={i,j,k}
        if any(mk<=d for d in done): continue
        (x1,y1),(x2,y2),(x3,y3)=P[i],P[j],P[k]
        d=2*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))
        if d==0:
            memb={t for t in range(n) if (P[t][0]-x1)*(y2-y1)-(P[t][1]-y1)*(x2-x1)==0}
            bl.append(('L',memb)); done.append(memb)
        else:
            a2=x1*x1+y1*y1;b2=x2*x2+y2*y2;c2=x3*x3+y3*y3
            cx=(a2*(y2-y3)+b2*(y3-y1)+c2*(y1-y2))/d; cy=(a2*(x3-x2)+b2*(x1-x3)+c2*(x2-x1))/d
            r=math.hypot(x1-cx,y1-cy)
            memb={t for t in range(n) if abs(math.hypot(P[t][0]-cx,P[t][1]-cy)-r)<1e-9}
            bl.append(('C',memb,cx,cy,r)); done.append(memb)
    return bl
P=[(0,2),(1,2),(2,2),(0,3),(2,3),(0,4),(1,4),(2,4)]
bl=blocks(P)
print(len(bl)); 
for b in bl: print(b)
