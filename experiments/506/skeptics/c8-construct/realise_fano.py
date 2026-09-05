"""Try to realise (numerically, random restarts) the structures that beat 17 combinatorially but
are excluded only by Sylvester-Gallai (a point in seven 4-blocks = derived Fano plane), plus the
Fano plane itself as 7 collinear-triple lines, and the Moebius-Kantor 8_3 as 8 lines.
If any of these were realisable, the lower bound would be refuted."""
import numpy as np, itertools
from realise import realise
sqs = [[0,1,2,3],[0,1,4,5],[2,3,4,5],[1,2,4,6],[0,3,4,6],[0,2,5,6],[1,3,5,6],[0,2,4,7],[1,3,4,7],[1,2,5,7],[0,3,5,7],[0,1,6,7],[2,3,6,7],[4,5,6,7]]
tests = {
  "Fano plane as 7 lines on 7 points": (7, [[0,1,3],[1,2,4],[2,3,5],[3,4,6],[4,5,0],[5,6,1],[6,0,2]], [[0,1,3],[1,2,4],[2,3,5],[3,4,6],[4,5,0],[5,6,1],[6,0,2]]),
  "Moebius-Kantor 8_3 as 8 lines": (8, [[i,(i+1)%8,(i+3)%8] for i in range(8)], [[i,(i+1)%8,(i+3)%8] for i in range(8)]),
  "SQS(8) (14 four-blocks) + 2 parallel lines": (8, sqs, [[0,1,2,3],[4,5,6,7]]),
  "SQS(8) minus one block, + 4 lines?": (8, sqs[1:], [[0,1,4,5],[2,3,6,7],[0,2,6],[1,3,7]]),
  "12 blocks (SQS minus 2 blocks sharing 2 pts) + 4 lines": (8, [b for b in sqs if b not in ([0,1,2,3],[0,1,4,5])], [[2,3,4,5],[0,1,6,7],[0,2,4],[1,3,5]]),
  "7 four-blocks through point 0 (derived Fano) alone": (8, [[0]+[a+1 for a in l] for l in [[0,1,3],[1,2,4],[2,3,5],[3,4,6],[4,5,0],[5,6,1],[6,0,2]]], []),
}
for name, (n, blocks, lines) in tests.items():
    print(name)
    P, nc, ncol = realise(n, blocks, lines, tries=300, seed=1)
    if P is not None:
        print("  REALISED?!", P, nc, ncol)
