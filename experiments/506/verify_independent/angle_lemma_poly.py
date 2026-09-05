import sympy as sp
m2,m3,m4,c2,c3,c4=sp.symbols('m2 m3 m4 c2 c3 c4')
L={1:(sp.Integer(0),sp.Integer(0)),2:(m2,c2),3:(m3,c3),4:(m4,c4)}
def P(i,j):
    mi,ci=L[i]; mj,cj=L[j]; d=mi-mj
    # homogeneous: (X, Y, W) with x=X/W, y=Y/W, W=d
    X=cj-ci; W=d; Y=mi*X+ci*W
    return (sp.expand(X),sp.expand(Y),sp.expand(W))
def concyc(pts):
    rows=[]
    for X,Y,W in pts:
        rows.append([sp.expand(X*X+Y*Y), sp.expand(X*W), sp.expand(Y*W), sp.expand(W*W)])
    return sp.Matrix(rows).det(method='berkowitz')
Q={1:[P(1,2),P(3,4),P(1,3),P(2,4)],2:[P(1,2),P(3,4),P(1,4),P(2,3)],3:[P(1,3),P(2,4),P(1,4),P(2,3)]}
A={}
for i,q in Q.items():
    d=sp.expand(concyc(q))
    f=sp.factor_list(d)
    print(f"det{i}: const {f[0]}")
    for g,e in f[1]: print("   ",g,"^",e)
    A[i]=[g for g,e in f[1] if g.has(m2) and g.has(m3) and g.has(m4) and not g.has(c2)]
    print("   essential:",A[i], flush=True)
print("A1-A2:",sp.expand(A[1][0]-A[2][0]))
print("A1-A3:",sp.expand(A[1][0]-A[3][0]))
