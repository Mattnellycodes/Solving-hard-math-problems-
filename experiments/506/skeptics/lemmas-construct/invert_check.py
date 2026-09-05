from fractions import Fraction as Fr
from count import analyse, blocks_of
import sys, ast
def invert(P, q):
    qx, qy = Fr(q[0]), Fr(q[1]); out = []
    for x, y in P:
        x, y = Fr(x), Fr(y)
        if (x, y) == (qx, qy): continue
        dx, dy = x - qx, y - qy; r2 = dx*dx + dy*dy
        out.append((qx + dx / r2, qy + dy / r2))
    return out
Q = ast.literal_eval(sys.argv[1])
Q = [(Fr(a), Fr(b)) for a, b in Q]
bl = blocks_of(Q)
print('Q blocks:', {k: sorted(v) for k, v in bl.items()})
for q in Q:
    P = invert(Q, q)
    print('invert about', q, '->', analyse(P))
