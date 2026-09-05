"""Targeted numerical realisation attacks on the structures the report declares unrealisable."""
import json, sys, itertools
from realise import attempt
tests = {}
# 1. n=6, 7 circles: complete quadrilateral (4 lines) + 3 concyclic opposite-pair complements
tests['n6_quadrilateral_7circles'] = dict(n=6, blocks=[(2,3,4,5),(0,1,4,5),(0,1,2,3)], lines=[(0,2,4),(0,3,5),(1,2,5),(1,3,4)])
# 1b. same but only two of the three concyclic (should be realisable: 8 circles) -- control
tests['n6_control_two_concyclic'] = dict(n=6, blocks=[(2,3,4,5),(0,1,4,5)], lines=[(0,2,4),(0,3,5),(1,2,5),(1,3,4)])
# 2. n=7: Fano complements as 4-blocks (no lines needed to beat 11? count = 35-21-l; need l>=4 for <=10; test blocks alone)
fano = [(0,1,2),(0,3,4),(0,5,6),(1,3,5),(1,4,6),(2,3,6),(2,4,5)]
comp = [tuple(sorted(set(range(7)) - set(l))) for l in fano]
tests['n7_fano_complement_4blocks'] = dict(n=7, blocks=comp, lines=[])
# 2b. control: 6 of the 7 complements
tests['n7_control_six_complements'] = dict(n=7, blocks=comp[:6], lines=[])
# 3. Fano plane as 3-point lines on 7 points (n=7, 7 lines) - classical impossibility
tests['n7_fano_lines'] = dict(n=7, blocks=[], lines=fano)
# 4. Mobius-Kantor (8_3) as lines on 8 points
mk = [tuple(sorted({i, (i+1)%8, (i+3)%8})) for i in range(8)]
tests['n8_MK_lines'] = dict(n=8, blocks=[], lines=mk)
# 4b. control: 7 of the 8 MK lines
tests['n8_control_MK_minus_line'] = dict(n=8, blocks=[], lines=mk[:7])
# 5. n=8: derived Fano at point 7: 7 four-blocks {7}+fano line (equivalent to 3 but as circles through a point)
tests['n8_fano_through_point'] = dict(n=8, blocks=[tuple(sorted(l + (7,))) for l in fano], lines=[])
# 6. n=8 cube structure with 4 lines? (report: l_max=3).  Cube: vertices = 3-bit strings; faces + diagonal planes.
V = list(itertools.product((0,1), repeat=3)); idx = {v: i for i, v in enumerate(V)}
faces = []
for ax in range(3):
    for val in (0,1):
        faces.append(tuple(sorted(idx[v] for v in V if v[ax] == val)))
diag = []
for ax in range(3):
    o = [a for a in range(3) if a != ax]
    for s in (0,1):
        diag.append(tuple(sorted(idx[v] for v in V if (v[o[0]] + v[o[1]] + s) % 2 == 0)))
cube = faces + diag
corner = [tuple(sorted(idx[tuple(v[j] ^ (j == a) for j in range(3))] for a in range(3))) for v in V]
tests['n8_cube_check_3lines'] = dict(n=8, blocks=cube, lines=[faces[0], corner[idx[(1,0,1)]], corner[idx[(0,1,1)]]])  # face z=0? plus two corner triples off it
if __name__ == '__main__':
    which = [a for a in sys.argv[1:] if not a.isdigit()] or list(tests)
    for name in which:
        t = tests[name]
        print('==', name, 'blocks', len(t['blocks']), 'lines', len(t['lines']), flush=True)
        r = attempt(t['n'], t['blocks'], t['lines'], tries=int(sys.argv[-1]) if sys.argv[-1].isdigit() else 60)
        print('   ->', 'REALISED' if r['success'] else 'not realised', 'residual', r['residual'], 'circles', r.get('circles'), flush=True)
        if r['success']:
            json.dump(r, open(f'realised_{name}.json', 'w'))
