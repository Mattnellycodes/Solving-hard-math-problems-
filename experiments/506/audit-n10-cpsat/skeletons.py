"""Auditor's own enumeration of 'skeletons': families of blocks of size >= 5 on 10 points, pairwise
sharing <= 2 points, satisfying the per-point Sylvester-Gallai cap sum_{B∋p} C(|B|-1,2) <= 35.
Up to isomorphism (S_10)."""
import itertools, json
from math import comb
from common import ClassCollector, N

big = [frozenset(c) for k in range(5, 10) for c in itertools.combinations(range(N), k)]


def ok_add(fam, B):
    for A in fam:
        if len(A & B) > 2:
            return False
    for p in B:
        s = sum(comb(len(A) - 1, 2) for A in fam if p in A) + comb(len(B) - 1, 2)
        if s > 35:
            return False
    return True


level = {0: [[]]}
allclasses = []
cur = ClassCollector()
cur.add([], [])
k = 0
while True:
    nxt = ClassCollector()
    for F, _, _ in cur.classes():
        for B in big:
            if B in F:
                continue
            if ok_add(F, B):
                nxt.add(sorted(F + [B], key=lambda s: (len(s), sorted(s))), [])
    if not nxt.reps:
        break
    k += 1
    cl = nxt.classes()
    print('level', k, ':', len(cl), 'classes')
    for F, _, _ in cl:
        allclasses.append([sorted(B) for B in F])
    cur = nxt
print('total skeletons (nonempty):', len(allclasses))
for F in allclasses:
    print(sorted(len(B) for B in F), F)
json.dump(allclasses, open('skeletons_audit.json', 'w'))
