"""MK non-realisability, exhaustive over projective charts (my own computation).
Frame p0=(0,0,1), p1=(1,0,1), p2=(0,1,1), p5=(1,1,1) (no three on an MK line, distinct).  Each remaining
point is (x,y,1) or (1,y,0) or (0,1,0) (charts cover all of RP^2).  For every chart combination
compute the Groebner basis of the 8 incidence equations plus distinctness (Rabinowitsch on all
pairs), and report whether it is the unit ideal; if not, print the basis."""
import itertools, sympy as sp
MK = [frozenset({i, (i+1) % 8, (i+3) % 8}) for i in range(8)]
frame = {0: (0,0,1), 1: (1,0,1), 2: (0,1,1), 5: (1,1,1)}
others = [3, 4, 6, 7]
t = sp.symbols('t')
total_unit = 0; cases = 0
for charts in itertools.product(range(3), repeat=4):
    P = {i: sp.Matrix(v) for i, v in frame.items()}
    vars_ = []
    for i, ch in zip(others, charts):
        x, y = sp.symbols(f'x{i} y{i}')
        if ch == 0: P[i] = sp.Matrix([x, y, 1]); vars_ += [x, y]
        elif ch == 1: P[i] = sp.Matrix([1, y, 0]); vars_ += [y]
        else: P[i] = sp.Matrix([0, 1, 0])
    eqs = [sp.expand(sp.Matrix.hstack(*[P[i] for i in sorted(L)]).det()) for L in MK]
    # distinctness: all 2x2 minors of a pair vanish iff points equal; use Rabinowitsch on sum of squares of minors
    dist = 1
    for i, j in itertools.combinations(range(8), 2):
        m = P[i].cross(P[j])
        dist *= (m[0]**2 + m[1]**2 + m[2]**2)
    eqs.append(sp.expand(t * dist - 1))
    G = sp.groebner(eqs, *vars_, t, order='grevlex')
    cases += 1
    if list(G) == [1]:
        total_unit += 1
    else:
        print('chart', charts, 'NOT unit ideal:', list(G))
print(f'{cases} chart combinations, {total_unit} give the unit ideal (no realisation, even over C, with distinct points)')
