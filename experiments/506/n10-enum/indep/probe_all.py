import json, math, time, sys
import numpy as np
import realise_numeric as rn
import mo
sg = json.load(open('n10_sg_t32.json'))['candidates']
rig = json.load(open('c2_rigidity.json')); ls = json.load(open('c2_linesets.json'))
STARTS = int(sys.argv[1]) if len(sys.argv) > 1 else 60
def run(name, blocks, lines, starts=STARTS):
    t0 = time.time()
    (r, P), nc, nl = rn.probe(10, blocks, lines, starts=starts, seed=1)
    print(f"{name}: {nc} concyclicity + {nl} collinearity equations; best residual over {starts} starts = {r:.3e}  [{time.time()-t0:.0f}s]", flush=True)
    return r, P
# control 1: the 17-circle 8-point configuration (cube structure with 3 lines) is realisable
cube = [[0,1,2,3],[0,1,4,5],[2,3,4,5],[0,2,4,6],[1,3,4,6],[1,2,5,6],[0,3,4,7],[0,2,5,7],[1,3,5,7],[0,1,6,7],[2,3,6,7],[4,5,6,7]]
(r, P), nc, nl = rn.probe(8, cube, [[0,1,2,3],[0,2,5,7],[1,3,4,6]], starts=30, seed=1)
print(f"control (n=8 cube + 3 lines, realisable): best residual {r:.3e}")
# control 2: C-2 alone (realisable: two regular pentagons)
r, P = run("C-2 alone", rig['blocks'], [])
if P is not None and r < 1e-12:
    print("   found realisation, points:", [(round(float(x), 4), round(float(y), 4)) for x, y in P])
# C-2 with each line-set orbit representative
for k, orb in enumerate(ls['orbits']):
    run(f"C-2 + line set orbit {k} (|L|={len(orb['lines'])})", rig['blocks'], orb['lines'])
# Miquel-excluded candidates 0 and 3 (blocks only)
for ci in (0, 3):
    run(f"candidate {ci} (Miquel-excluded) blocks only", sg[ci]['blocks'], [])
