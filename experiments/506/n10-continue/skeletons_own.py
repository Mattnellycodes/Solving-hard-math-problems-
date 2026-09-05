"""Own enumeration of the big-block skeletons for n = 10: isomorphism classes of families of blocks of size
5..9 on 10 points, pairwise sharing <= 2 points, subject only to the per-point Sylvester–Gallai cap
sum_{B ∋ p} C(|B|-1, 2) <= C(9,2) - 1 = 35 (derived 9-point set has >= 1 ordinary line).
BFS by number of blocks with isomorphism reduction (invariant + VF2).  Compared with n10-cpsat/v2/skeletons.json.
"""
import itertools, json, math, sys, time
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-continue')
from common import ClassStore, fs, N

CAP = math.comb(9, 2) - 1
cands = [fs(c) for k in range(5, 10) for c in itertools.combinations(range(N), k)]

def ok(fam, B):
    if any(len(B & C) > 2 for C in fam):
        return False
    for p in range(N):
        if sum(math.comb(len(C) - 1, 2) for C in fam + [B] if p in C) > CAP:
            return False
    return True

t0 = time.time()
allc = []
frontier = [[]]
level = 0
while frontier:
    level += 1
    store = ClassStore(N)
    for fam in frontier:
        for B in cands:
            if B in fam or not ok(fam, B):
                continue
            store.add(fam + [B], [])
    frontier = [c['F'] for c in store.classes]
    print(f"level {level}: {len(frontier)} classes  [{time.time()-t0:.0f}s]", flush=True)
    allc += frontier
print("total classes:", len(allc))
recs = []
for fam in allc:
    recs.append({'sizes': sorted(len(B) for B in fam), 'blocks': sorted(sorted(B) for B in fam)})
    print(recs[-1]['sizes'], recs[-1]['blocks'])
json.dump(recs, open('skeletons_own.json', 'w'), indent=0)
# compare with v2's 18 skeletons (cases A, B) + the big-block ones
v2 = json.load(open('/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2/skeletons.json'))
store = ClassStore(N)
for r in v2:
    store.add([fs(B) for B in r['blocks']], [])
n_v2 = len(store.classes)
matched = 0
extra = []
for fam in allc:
    idx, new = store.add(fam, [])
    if new:
        extra.append(sorted(len(B) for B in fam))
    else:
        matched += 1
print(f"v2 skeleton classes: {n_v2}; own classes matching a v2 class: {matched}; own classes not in v2: {extra}")
print(f"done [{time.time()-t0:.0f}s]")
