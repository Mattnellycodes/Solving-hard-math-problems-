"""Enumerate, up to isomorphism, all families of 5-subsets of {0..9} that contain {0,1,2,3,4},
pairwise share <= 2 points and have every point in <= 3 members (proved bound, check_local_bounds.py).
These are the possible sets of 5-blocks of a case-B structure (largest block 5).
Output: skeletons.json = list of {'k': size, 'blocks': [...]}.
"""
import itertools, json, time, sys
from collections import Counter
import networkx as nx

N = 10
PTS = range(N)
FIVES = [frozenset(c) for c in itertools.combinations(PTS, 5)]
BASE = frozenset(range(5))


def graph(F):
    G = nx.Graph()
    for p in PTS:
        G.add_node(('p', p), c=0)
    for i, B in enumerate(F):
        G.add_node(('b', i), c=1)
        for p in B:
            G.add_edge(('b', i), ('p', p))
    return G


def inv(F):
    deg = Counter(p for B in F for p in B)
    pair = Counter()
    for B1, B2 in itertools.combinations(F, 2):
        pair[len(B1 & B2)] += 1
    return (len(F), tuple(sorted(deg[p] for p in PTS)), tuple(sorted(pair.items())),
            tuple(sorted(tuple(sorted(deg[q] for q in B)) for B in F)),
            nx.weisfeiler_lehman_graph_hash(graph(F), node_attr='c', iterations=3))


nm = nx.algorithms.isomorphism.categorical_node_match('c', 0)

def main(kmax=6):
    t0 = time.time()
    levels = {1: [[BASE]]}
    allrecs = [{'k': 1, 'blocks': [sorted(BASE)]}]
    for k in range(2, kmax + 1):
        buckets = {}
        reps = []
        cand_count = 0
        for F in levels[k - 1]:
            deg = Counter(p for B in F for p in B)
            for B in FIVES:
                if B in F or any(len(B & C) > 2 for C in F) or any(deg[p] >= 3 for p in B):
                    continue
                F2 = F + [B]
                cand_count += 1
                key = inv(F2)
                G = graph(F2)
                found = False
                for idx in buckets.get(key, []):
                    if nx.is_isomorphic(reps[idx][1], G, node_match=nm):
                        found = True; break
                if not found:
                    buckets.setdefault(key, []).append(len(reps))
                    reps.append((F2, G))
        levels[k] = [r[0] for r in reps]
        allrecs += [{'k': k, 'blocks': sorted(sorted(B) for B in F)} for F in levels[k]]
        print(f"k={k}: {len(levels[k])} classes ({cand_count} candidates, {time.time()-t0:.0f}s)", flush=True)
    json.dump(allrecs, open('skeletons.json', 'w'))
    print("total", len(allrecs))

if __name__ == '__main__':
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
