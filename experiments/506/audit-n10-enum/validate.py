"""V1: orderly generation (stage1 + stage2, incl. empty big family with G=S_n) reproduces the brute-force
isomorphism classes of ALL labelled families of rich blocks (sizes 4..n-1, pairwise <=2, derived cap) for n=6,7.
V2: stage-2 canonicity test agrees with brute force on random families (n=7,8)."""
import sys, itertools, random, time
import numpy as np
from math import comb
from mob import *

def brute_classes(n, capD):
    subs = [mask_of(c) for k in range(4, n) for c in itertools.combinations(range(n), k)]
    subs.sort()
    fams = []
    def rec(F, cov, start):
        fams.append(list(F))
        for i in range(start, len(subs)):
            m = subs[i]
            if any(popcount(m & b) > 2 for b in F): continue
            k = popcount(m); add = comb(k-1, 2)
            if any(cov[p] + add > capD for p in bits(m)): continue
            cov2 = cov[:]
            for p in bits(m): cov2[p] += add
            rec(F + [m], cov2, i + 1)
    rec([], [0]*n, 0)
    return fams

for n in (6, 7):
    capD, capL = caps(n, 'sg')
    PT = PermTable(n)
    t0 = time.time()
    fams = brute_classes(n, capD)
    classes = set()
    for F in fams:
        if F:
            classes.add(canon_full(PT, [bits(m) for m in F]))
    print(f"n={n}: {len(fams)-1} labelled nonempty families, {len(classes)} brute-force classes [{time.time()-t0:.1f}s]")
    # orderly
    got = set()
    seen_list = []
    for F, aut in stage1(n, capD, range(5, n), PT=PT, verbose=False):
        S = Stage2(n, F, aut, capD, capL, 0)
        S.collect_only = True
        S.run()
        for fam in S.candidates:
            if fam:
                seen_list.append(tuple(fam))
    print(f"   orderly visited {len(seen_list)} nonempty families; distinct={len(set(seen_list))}")
    got = set(canon_full(PT, [bits(m) for m in fam]) for fam in seen_list)
    ok = (got == classes) and len(seen_list) == len(set(seen_list)) == len(classes)
    print(f"   canonical forms match brute force: {got == classes}; each class exactly once: {len(seen_list)==len(classes)};  V1 {'OK' if ok else 'FAIL'}")

# V2: canonicity test vs brute force
random.seed(1)
for n in (7, 8):
    PT = PermTable(n)
    capD, capL = caps(n, 'sg')
    bad = 0; tests = 0
    for trial in range(300):
        # random big family: 0 or 1 or 2 blocks of size 5
        big = []
        for _ in range(random.choice([0, 1, 1, 2])):
            m = mask_of(random.sample(range(n), 5))
            if all(popcount(m & b) <= 2 for b in big): big.append(m)
        big.sort()
        # brute-force canonicalise big family so that stage-2 canonical form is well defined
        imgs = PT.image_masks(big) if big else None
        if big:
            order = np.lexsort(imgs.T[::-1]); big = [int(x) for x in imgs[order[0]]]
            less, equal = rows_less_or_equal(PT.image_masks(big), big); G = PT.P[equal]
        else:
            G = PT.P
        S = Stage2(n, big, G, capD, capL, 0)
        # random compatible 4-family
        idx = list(range(S.K)); random.shuffle(idx); F = []
        for i in idx:
            if all(popcount(S.cmask[i] & S.cmask[j]) <= 2 for j in F): F.append(i)
            if len(F) >= random.randint(1, 8): break
        F.sort()
        mine = S.canonical(F)
        # brute force: min over G of sorted code list
        codes = sorted(S.ccode[i] for i in F)
        pts = [S.cands[i] for i in F]
        best = None
        for g in G:
            im = sorted(code4(tuple(sorted(int(g[p]) for p in c))) for c in pts)
            if best is None or im < best: best = im
        truth = (best == codes)
        tests += 1
        if truth != mine: bad += 1
    print(f"V2 n={n}: {tests} canonicity tests, {bad} disagreements  {'OK' if bad==0 else 'FAIL'}")
