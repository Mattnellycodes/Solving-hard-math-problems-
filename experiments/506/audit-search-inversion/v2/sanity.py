"""Sanity checks of count506 on known certificates (theory/REPORT.md) and vs the repo counter on random sets."""
import sys, random
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506')
from count506 import analyse, f
from math import comb
known = {6: ([(0,0),(20,0),(5,15),(5,0),(2,6),(10,10)], 8),
         7: ([(0,0),(20,0),(5,15),(5,0),(2,6),(10,10),(5,5)], 11),
         8: ([(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)], 17)}
for n, (pts, want) in known.items():
    r = analyse(pts)
    print(f"known n={n}: my count={r['circles']} expected={want} lines={r['line_sizes']} circle_sizes={r['circle_sizes']} OK={r['circles']==want}")
# cross-check vs repo exact counter on random small integer sets (consistency only)
try:
    import circles_exact as ce
    fn = [getattr(ce, a) for a in dir(ce) if 'count' in a.lower() and callable(getattr(ce, a))]
    print("repo counter functions:", [g.__name__ for g in fn])
except Exception as e:
    print("repo counter import failed:", e)
random.seed(1)
mism = 0
for t in range(300):
    n = random.randint(4, 9)
    pts = set()
    while len(pts) < n:
        pts.add((random.randint(-3, 3), random.randint(-3, 3)))
    pts = sorted(pts)
    mine = analyse(pts)['circles']
    theirs = ce.count_circles(pts) if hasattr(ce, 'count_circles') else None
    if isinstance(theirs, (tuple, list, dict)):
        theirs = theirs[0] if isinstance(theirs, (tuple, list)) else theirs.get('circles')
    if theirs != mine:
        mism += 1
        if mism < 5: print("MISMATCH", pts, mine, theirs)
print("random cross-check mismatches:", mism)
