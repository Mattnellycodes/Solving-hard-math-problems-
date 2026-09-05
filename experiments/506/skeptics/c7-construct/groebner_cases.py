import sympy as sp
mb,mc,md,t=sp.symbols('mb mc md t')
Phi = mb*mc*t - mb*md*t + mb*md - mc*md
# 1. Phi = 0  <=>  L_b, L_c, L_d concurrent (P4 = P6)?
P4=(mb*t/(mb-mc), mb*mc*t/(mb-mc)); P6=(-md/(mc-md), -mc*md/(mc-md))
print("P4==P6 numerator:", sp.factor(sp.numer(sp.together(P4[0]-P6[0]))), "|", sp.factor(sp.numer(sp.together(P4[1]-P6[1]))))
# 2. component 2 from solve: does it satisfy Phi=0?
comp2={mb:(md-mc)/(mc*md+1), t: md*(-mc**2*md-2*mc+md)/(mc**2-2*mc*md+md**2)}
print("Phi on component 2:", sp.simplify(Phi.subs(comp2)))
# 3. vertical-line cases (slope parametrisation misses them). L_a: y=0 always (rotation used).
def run_case(name, lines):
    # lines: dict name->(kind, params): ('slope',m,c) y=m x+c  or ('vert',x0)
    def inter(L1,L2):
        if L1[0]=='vert' and L2[0]=='vert': return None
        if L1[0]=='vert': x=L1[1]; return (x, L2[1]*x+L2[2])
        if L2[0]=='vert': x=L2[1]; return (x, L1[1]*x+L1[2])
        x=(L2[2]-L1[2])/(L1[1]-L2[1]); return (sp.simplify(x), sp.simplify(L1[1]*x+L1[2]))
    La,Lb,Lc,Ld=lines['a'],lines['b'],lines['c'],lines['d']
    P={0:inter(La,Lc),1:inter(La,Ld),3:inter(La,Lb),4:inter(Lb,Lc),5:inter(Lb,Ld),6:inter(Lc,Ld)}
    def cyc(ids):
        M=sp.Matrix([[P[i][0],P[i][1],P[i][0]**2+P[i][1]**2,1] for i in ids]); return sp.factor(sp.numer(sp.together(M.det())))
    Q=[cyc([0,1,4,5]),cyc([1,3,4,6]),cyc([0,3,5,6])]
    print("CASE",name); 
    for q in Q: print("   ",q)
    vars_=sorted(set().union(*[q.free_symbols for q in Q]),key=str)
    print("   solve:",sp.solve(Q,vars_,dict=True))
    print("   points:",P)
# L_c vertical (x=0), P0=(0,0), P1=(1,0)
run_case("Lc vertical", {'a':('slope',0,0),'b':('slope',mb,-mb*t),'c':('vert',0),'d':('slope',md,-md)})
run_case("Ld vertical", {'a':('slope',0,0),'b':('slope',mb,-mb*t),'c':('slope',mc,0),'d':('vert',1)})
run_case("Lb vertical", {'a':('slope',0,0),'b':('vert',t),'c':('slope',mc,0),'d':('slope',md,-md)})
