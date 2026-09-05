import sys, re
from invert_attack import best_inversion
best = 10**9; seen=0
for line in open(sys.argv[1]):
    if not line.startswith("SET"): continue
    pts = [(int(a), int(b)) for a, b in re.findall(r"\((-?\d+),(-?\d+)\)", line)]
    r = best_inversion(pts); seen+=1
    if r is None: continue
    nb, nl, deg, O, cnt = r
    if cnt < best: best = cnt; print("new best", cnt, "from", pts, "|B|=", nb, "deg", deg, "O", O)
print("sets processed", seen, "best count after optimal inversion:", best)
