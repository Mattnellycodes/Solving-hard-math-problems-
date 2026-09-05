"""The 1-parameter family of 8-point sets with 17 circles: A=(0,3), B=(b,0), C=(3/b,0), 0<b<sqrt3.
H=(0,-1) is the orthocentre (bc = 3 = a*h), O=(0,1) = midpoint of AH = reflection of H in BC lies on
the altitude, the circumcircle and the nine-point circle.  Invert {A,B,C,H,D,E,F} in O and count."""
from fractions import Fraction as Fr
from certify import analyse, f
def foot(P, Q, R):
    """foot of the perpendicular from P onto line QR (exact)."""
    dx, dy = R[0]-Q[0], R[1]-Q[1]
    t = ((P[0]-Q[0])*dx + (P[1]-Q[1])*dy) / (dx*dx + dy*dy)
    return (Q[0] + t*dx, Q[1] + t*dy)
def config(b):
    b = Fr(b); A = (Fr(0), Fr(3)); B = (b, Fr(0)); C = (3/b, Fr(0)); H = (Fr(0), Fr(-1)); O = (Fr(0), Fr(1))
    D = foot(C, A, B); E = foot(B, A, C); F = foot(A, B, C)
    assert F == (0, 0)
    def inv(p):
        dx, dy = p[0]-O[0], p[1]-O[1]; d2 = dx*dx + dy*dy
        return (O[0] + dx/d2, O[1] + dy/d2)
    return [inv(p) for p in (A, B, C, H, D, E, F)] + [O]
for b in [Fr(1), Fr(1,2), Fr(3,2), Fr(2,3), Fr(5,4), Fr(1,3), Fr(7,5)]:
    P = config(b)
    print(f"b={b}: n=8 -> circles, 3+lines, circle sizes, line sizes, degenerate = {analyse(P)}")
