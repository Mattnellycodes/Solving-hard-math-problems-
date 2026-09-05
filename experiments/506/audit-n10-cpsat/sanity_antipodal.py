"""Sanity test of the filters: realisable configurations must pass every filter (own exact structure
extraction with Fractions).  Includes the 33-circle antipodal configuration for n = 10."""
from fractions import Fraction as Fr
import itertools
from common import run_filters, count_of, check_structure, D_of

def structure(pts):
    n = len(pts)
    blocks = {}
    for i, j, k in itertools.combinations(range(n), 3):
        (x1,y1),(x2,y2),(x3,y3) = pts[i],pts[j],pts[k]
        det = (x2-x1)*(y3-y1) - (y2-y1)*(x3-x1)
        if det == 0:
            a = y2-y1; b = x1-x2; c = -(a*x1+b*y1)
            key = ('L',) + tuple(t/ (a if a else b) for t in (a,b,c))
        else:
            A = x2-x1; B = y2-y1; C = x3-x1; Dd = y3-y1
            E = A*(x1+x2)+B*(y1+y2); F = C*(x1+x3)+Dd*(y1+y3)
            G = 2*(A*(y3-y2)-B*(x3-x2))
            cx = (Dd*E-B*F)/G; cy = (A*F-C*E)/G
            key = ('C', cx, cy, (cx-x1)**2+(cy-y1)**2)
        blocks.setdefault(key, set()).update((i,j,k))
    F = [sorted(v) for k, v in blocks.items() if len(v) >= 4]
    L = [sorted(v) for k, v in blocks.items() if k[0] == 'L']
    circles = sum(1 for k, v in blocks.items() if k[0] == 'C')
    return F, L, circles

base = [(Fr(1),Fr(0)),(Fr(-1),Fr(0)),(Fr(0),Fr(1)),(Fr(0),Fr(-1)),
       (Fr(3,5),Fr(4,5)),(Fr(-3,5),Fr(-4,5)),(Fr(3,5),Fr(-4,5)),(Fr(-3,5),Fr(4,5)),(Fr(0),Fr(0))]
# 9 concyclic points (4 antipodal pairs + one more) + the centre: Lemma A extremal, 33 circles
extra = (Fr(5,13), Fr(12,13))
F, L, circ = structure(base + [extra])
check_structure(F, L)
print('best extra', extra, 'circles', circ, 'formula', count_of(F, L), 'D', D_of(F), 'l', len(L))
print('F =', F); print('L =', L)
print('  violations (proved filters):', run_filters(F, L))
print('  violations (+cited):', run_filters(F, L, cited=True))
# a second realisable set: 3x3 grid + centre-ish point
grid = [(Fr(x),Fr(y)) for x in range(3) for y in range(3)] + [(Fr(5),Fr(7))]
F, L, circ = structure(grid); check_structure(F, L)
print('grid+1: circles', circ, 'violations', run_filters(F, L), run_filters(F, L, cited=True))
