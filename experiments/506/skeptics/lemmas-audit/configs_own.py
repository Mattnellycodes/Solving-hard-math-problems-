"""lemmas-audit: independent combinatorial facts used by the report.
 (a) every family of 8 triples on 8 points pairwise sharing <= 1 point is (8_3) = Moebius-Kantor;
     more generally every family of 8 sets of size >= 3 on 8 points pairwise sharing <= 1 point.
 (b) 7 triples on 7 points pairwise sharing <= 1 point = Fano (unique); 8 impossible.
 (c) 6 triples on 7 points pairwise sharing <= 1 point = Fano minus a line (unique).
 (d) max number of pairwise <=1-sharing triples on 6 points is 4, on 5 points 2, on 9 points 12 (STS(9)).
 (e) 12 four-subsets of an 8-set pairwise sharing <= 2 points: classify (expect: exactly one class,
     AG(3,2) minus a parallel class = 'cube').  Also 7 four-subsets of a 7-set pairwise sharing <= 2
     = Fano complements (unique).
Everything by brute force / canonical forms computed with all permutations.
"""
import itertools
from math import comb

def canon(sets, n):
    best = None
    for s in itertools.permutations(range(n)):
        img = tuple(sorted(tuple(sorted(s[i] for i in t)) for t in sets))
        if best is None or img < best:
            best = img
    return best

def enum_partial_linear(n, k, size, maxshare=1):
    """all families of `size` k-subsets of [n] pairwise sharing <= maxshare points, up to iso."""
    subs = [frozenset(c) for c in itertools.combinations(range(n), k)]
    classes = set()
    count_labelled = 0
    def rec(F, start):
        nonlocal count_labelled
        if len(F) == size:
            count_labelled += 1
            classes.add(canon(F, n)); return
        for j in range(start, len(subs)):
            t = subs[j]
            if all(len(t & u) <= maxshare for u in F):
                rec(F + [t], j + 1)
    rec([], 0)
    return classes, count_labelled

def max_packing(n, k, maxshare=1):
    subs = [frozenset(c) for c in itertools.combinations(range(n), k)]
    best = 0
    def rec(F, start):
        nonlocal best
        best = max(best, len(F))
        for j in range(start, len(subs)):
            t = subs[j]
            if all(len(t & u) <= maxshare for u in F):
                rec(F + [t], j + 1)
    rec([], 0)
    return best

MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
FANO = [frozenset(t) for t in [{0,1,2},{0,3,4},{0,5,6},{1,3,5},{1,4,6},{2,3,6},{2,4,5}]]

if __name__ == "__main__":
    # (a) 8 triples on 8 points -- fix the first triple {0,1,2} WLOG is fine for iso classes, but do it fully:
    classes, lab = enum_partial_linear(8, 3, 8)
    print(f"(a) 8 triples on 8 points pairwise <=1: {lab} labelled families, {len(classes)} iso class(es); "
          f"MK in classes: {canon(MK, 8) in classes}; 8!/|Aut(MK)| = {40320 // lab if lab else None}")
    # also with mixed sizes >= 3: check impossibility of a 4-set among 8 pairwise <=1 sets on 8 points
    subs = [frozenset(c) for k in (3, 4, 5) for c in itertools.combinations(range(8), k)]
    found = 0
    def rec(F, start):
        global found
        if len(F) == 8:
            if any(len(t) >= 4 for t in F):
                found += 1
            return
        for j in range(start, len(subs)):
            t = subs[j]
            if all(len(t & u) <= 1 for u in F):
                rec(F + [t], j + 1)
    rec([], 0)
    print(f"    families of 8 pairwise<=1 sets of size>=3 on 8 points containing a set of size >=4: {found}")
    # (b)
    classes7, lab7 = enum_partial_linear(7, 3, 7)
    print(f"(b) 7 triples on 7 points: {lab7} labelled, {len(classes7)} class(es); Fano: {canon(FANO, 7) in classes7}; "
          f"8 triples possible: {max_packing(7, 3) >= 8} (max packing {max_packing(7, 3)})")
    # (c)
    classes6, lab6 = enum_partial_linear(7, 3, 6)
    print(f"(c) 6 triples on 7 points pairwise <=1: {len(classes6)} iso class(es) (expect 1 = Fano minus a line)")
    # (d)
    print(f"(d) max packings of triples: n=5: {max_packing(5,3)}, n=6: {max_packing(6,3)}, n=8: {max_packing(8,3)}")
    # (e) 7 four-sets of a 7-set pairwise sharing <= 2
    c74, l74 = enum_partial_linear(7, 4, 7, maxshare=2)
    print(f"(e) 7 four-subsets of a 7-set pairwise <=2: {len(c74)} class(es); 8 possible: {max_packing(7,4,2) >= 8}")
    # 12 four-subsets of an 8-set pairwise sharing <= 2 -- use degree bound (<= 6 per point is automatic:
    # 7 four-blocks through p would need 7 pairwise-disjoint... no: pairwise <=1 triples on 7 points -> max 7 (Fano),
    # so degree <= 7 combinatorially).  Enumerate with first block {0,1,2,3} fixed and canonical forms at leaves.
    subs4 = [frozenset(c) for c in itertools.combinations(range(8), 4)]
    leaves = set(); nlab = 0
    def rec12(F, start):
        global nlab
        if len(F) == 12:
            nlab += 1
            leaves.add(canon(F, 8)); return
        # potential
        if len(F) + (len(subs4) - start) < 12:
            return
        for j in range(start, len(subs4)):
            t = subs4[j]
            if all(len(t & u) <= 2 for u in F):
                rec12(F + [t], j + 1)
    rec12([subs4[0]], 1)
    print(f"(e) 12 four-subsets of an 8-set pairwise <=2 containing {sorted(subs4[0])}: {nlab} labelled, {len(leaves)} iso class(es)")
    for cl in leaves:
        degs = [sum(1 for t in cl if p in t) for p in range(8)]
        print("     class:", cl, "degrees", degs)
    print(f"    max packing of 4-sets on 8 points pairwise<=2: {max_packing(8,4,2)} (SQS(8) has 14)")
