"""(a) Every family of 8 triples on 8 points pairwise sharing <= 1 point is isomorphic to the
Möbius-Kantor configuration {i, i+1, i+3 mod 8}  (so t3(8) <= 7 follows from (8_3) non-real).
(b) Every family of 11 triples on 9 points pairwise sharing <= 1 point has an 8-subset carrying 8
of the triples (an (8_3)), so t3(9) <= 10.  Own orderly BFS (canonical dedupe per level)."""
from itertools import combinations
from lib import mask, popcount, canon, mk_violation

MK = [mask(((i) % 8, (i + 1) % 8, (i + 3) % 8)) for i in range(8)]
mk_cf = canon(MK, 8)[0]


def packing_classes(n, k):
    """Isomorphism classes (canonical forms) of families of k triples on n points, pairwise
    sharing <= 1 point, level by level."""
    tri = [mask(c) for c in combinations(range(n), 3)]
    level = {(): None}
    for step in range(k):
        nxt = {}
        for fam in level:
            for t in tri:
                if t in fam or any(popcount(t & b) >= 2 for b in fam):
                    continue
                nxt[canon(list(fam) + [t], n)[0]] = None
        level = nxt
        print(f"  n={n}: {step + 1} triples: {len(level)} classes", flush=True)
    return list(level)


P8 = packing_classes(8, 8)
print("classes of 8-triple packings on 8 points:", len(P8), "== MK:", set(P8) == {mk_cf})
assert set(P8) == {mk_cf}

P9 = packing_classes(9, 11)
print("classes of 11-triple packings on 9 points:", len(P9))
bad = [f for f in P9 if mk_violation(list(f), (1 << 9) - 1) is None]
print("classes WITHOUT an (8_3) sub-configuration:", len(bad))
assert not bad
print("t3(9) <= 10 follows from the (8_3) lemma")
