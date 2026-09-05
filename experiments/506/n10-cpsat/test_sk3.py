import itertools, sys
sys.path.insert(0, '.')
from realise import try_realise
from filters import full_report, bundle_violations
def build(bij):
    S = list(range(5)); R = list(range(5, 10))
    M = {v: [frozenset({S[(v-1)%5], S[(v+1)%5]}), frozenset({S[(v-2)%5], S[(v+2)%5]})] for v in range(5)}
    N = {w: [frozenset({R[(w-1)%5], R[(w+1)%5]}), frozenset({R[(w-2)%5], R[(w+2)%5]})] for w in range(5)}
    F = [frozenset(S), frozenset(R)]
    L = [frozenset(S)]
    for v in range(5):
        w = bij[v]
        for ab in M[v]:
            for rr in N[w]:
                F.append(ab | rr)
        for rr in N[w]:
            L.append(rr | {v})
    return F, L
for bij in [lambda v: v, lambda v: (2*v) % 5, lambda v: (-v) % 5]:
    F, L = build([bij(v) for v in range(5)])
    assert all(len(B1 & B2) <= 2 for B1, B2 in itertools.combinations(F, 2))
    assert all(len(S1 & S2) <= 1 for S1, S2 in itertools.combinations(L, 2))
    print("structure: blocks", len(F), "lines", len(L), "count", 120 - (9*2 + 3*20) - len(L))
    full_report(F, L)
    print("  bundle violations:", len(bundle_violations(F, L)))
    r = try_realise(F, L, tries=150, seed=1, verbose=True)
    print("  realised:", r['realised'], r.get('best_residual'))
