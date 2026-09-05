"""Compare my theorem-only survivor with the survivor actually stored by n10-continue
(runs/post_all.json, index 12) -- the block list printed in their REPORT/summary is mis-transcribed
(it contains the 3-line [2,4,6] inside the 4-block [2,4,6,7])."""
import json
from lib import mask, bits, canon, mobius_blocks, circles_count, run_filters, check_structure
mine = json.load(open("runs/survivors.json"))
theirs = json.load(open("../n10-continue/runs/post_all.json"))[12]
F2 = [mask(b) for b in theirs["blocks"]]; L2 = [mask(l) for l in theirs["lines"]]
check_structure(F2, L2)
print("their stored survivor (post_all #12): circles =", circles_count(F2, L2), "my filters on it:", run_filters(F2, L2))
cf2 = canon(mobius_blocks(F2, L2) + [1 << 10], 11)[0]
for s in mine:
    F = [mask(b) for b in s["F"]]; L = [mask(l) for l in s["L"]]
    cf = canon(mobius_blocks(F, L) + [1 << 10], 11)[0]
    print("my survivor", s["file"], s["index"], "isomorphic to theirs (inf fixed):", cf == cf2)
# the printed version in their REPORT
rep_F = [[0,1,2,3,4],[5,6,7,8,9],[0,1,6,9],[0,1,7,8],[0,2,5,9],[0,2,6,8],[0,3,5,7],[0,3,8,9],[0,4,5,6],[0,4,7,9],[1,2,5,8],[1,2,7,9],[1,3,5,9],[1,3,6,7],[1,4,5,7],[1,4,6,8],[2,3,5,6],[2,3,7,8],[2,4,6,7],[2,4,8,9],[3,4,5,8],[3,4,6,9]]
rep_L = [[0,1,5],[0,4,8],[0,6,7],[1,3,7],[1,8,9],[2,3,9],[2,4,6],[2,7,8],[3,5,6],[4,5,9]]
bad = [(l, b) for l in rep_L for b in rep_F if set(l) <= set(b)]
print("REPORT-printed survivor: 3-lines inside blocks:", bad)
