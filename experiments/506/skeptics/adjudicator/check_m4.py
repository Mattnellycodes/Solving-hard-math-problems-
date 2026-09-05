"""Check every family found by enum_big for n=9, m=4 (no 5-block): some point must lie in >= 8
four-blocks (then its derived structure is an (8_3), impossible over R).  Also histogram (D,ell,count)."""
import re, sys
from collections import Counter
fn = sys.argv[1]
hist = Counter(); bad = 0; tot = 0; mindeg = 99
for line in open(fn):
    if not line.startswith('FOUND'):
        continue
    tot += 1
    D = int(re.search(r'D=(\d+)', line).group(1)); e = int(re.search(r'ell=(\d+)', line).group(1)); c = int(re.search(r'count=(-?\d+)', line).group(1))
    blocks = [[int(x) for x in b.split(',')] for b in re.findall(r'\[([0-9,]+)\]', line)]
    deg = [sum(1 for B in blocks if p in B) for p in range(9)]
    md = max(deg); mindeg = min(mindeg, md)
    hist[(D, e, c, md)] += 1
    if md < 8:
        bad += 1
print(fn, "families:", tot, "histogram (D,ell,count,maxdeg):", dict(hist))
print("families with NO point in >= 8 four-blocks:", bad, "; minimum over families of the max degree:", mindeg)
print("END line:", [l.strip() for l in open(fn) if l.startswith('END')])
