"""Enumerate the 'skeletons' = isomorphism classes of families of big blocks (sizes 5 and 6) on 10 points that
can occur in a configuration with <= 32 circles whose largest block has size 5 or 6:
  case A: contains the 6-block {0..5}; other blocks of size 5 or 6 pairwise sharing <= 2 points;
  case B: contains the 5-block {0..4}, no 6-block; 5-blocks pairwise sharing <= 2 points.
Per-point bounds (facts.py F1, pure packing of the derived structure on 9 points): <= 2 six-blocks and
<= 3 five-blocks through a point.  Output skeletons.json.
"""
import itertools, json, sys
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2')
from model import ClassStore

N = 10
PTS = range(N)

def enumerate_case(case):
    if case == 'A':
        base = frozenset(range(6)); cands = [frozenset(c) for k in (5, 6) for c in itertools.combinations(PTS, k)]
    else:
        base = frozenset(range(5)); cands = [frozenset(c) for c in itertools.combinations(PTS, 5)]
    store = ClassStore(N)
    store.add([base], [])
    frontier = [[base]]
    level = 1
    while frontier:
        level += 1
        new_store = ClassStore(N)
        for fam in frontier:
            for B in cands:
                if B in fam or any(len(B & C) > 2 for C in fam):
                    continue
                fam2 = fam + [B]
                deg6 = {p: sum(1 for C in fam2 if p in C and len(C) == 6) for p in PTS}
                deg5 = {p: sum(1 for C in fam2 if p in C and len(C) == 5) for p in PTS}
                if any(deg6[p] > 2 or deg5[p] > 3 for p in PTS):
                    continue
                new_store.add(fam2, [])
        frontier = [c['F'] for c in new_store.classes]
        for c in new_store.classes:
            store.add(c['F'], [])
        print(f"case {case}: level {level}: {len(frontier)} classes", flush=True)
    return [{'case': case, 'k': len(c['F']), 'blocks': sorted(sorted(B) for B in c['F'])} for c in store.classes]

if __name__ == '__main__':
    recs = enumerate_case('A') + enumerate_case('B')
    for i, r in enumerate(recs):
        print(i, r['case'], r['k'], r['blocks'])
    json.dump(recs, open('skeletons.json', 'w'), indent=0)
    print("total skeletons:", len(recs))
