from realise import attempt
from count import blocks_of
for P in ([(0,0),(20,0),(5,15),(5,0),(2,6),(10,10)], [(0,0),(0,5),(0,10),(0,15),(3,6),(3,9),(5,5),(5,10)],
          [(0,0),(1,0),(-1,0),(0,1),(0,-1),(3,4),(-3,-4),(3,-4),(-3,4)]):
    bl = blocks_of(P)
    L = [tuple(sorted(s)) for k, s in bl.items() if k[0]=='L' and len(s)>=3]
    C = [tuple(sorted(s)) for k, s in bl.items() if k[0]=='C' and len(s)>=4]
    print('n=%d lines' % len(P), L, 'circles>=4', C)
    r = attempt(len(P), C, L, tries=100); print(r['success'], r.get('circles'))
