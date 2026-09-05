#!/usr/bin/env python3
"""Validation of the canonical-form machinery of oe.py against brute force.
V1: for n = 5, 6, 7 (modes none/table) the set of families visited by the two-phase orderly
    generator equals the set of brute-force canonical forms of ALL labelled families of rich
    blocks (sizes 4..n-1, pairwise <= 2 common points, cap C2).
V2: is_canonical_marked agrees with canon_bruteforce on random families (n = 7, 8).
"""
import sys, time, itertools, random
from math import comb
import numpy as np
from oe import *

def labelled_families(n, caps):
    blocks = []
    for s in range(4, n):
        blocks += [list_to_mask(c) for c in itertools.combinations(range(n), s)]
    blocks.sort()
    K = len(blocks)
    fams = []
    def rec(start, fam, used):
        fams.append(tuple(fam))
        for i in range(start, K):
            b = blocks[i]
            if any(popcount(b & B) > 2 for B in fam):
                continue
            w = comb(popcount(b) - 1, 2)
            pts = mask_to_list(b, n)
            if any(used[p] + w > caps.capD for p in pts):
                continue
            if popcount(b) == 4 and any(sum(1 for B in fam if popcount(B) == 4 and B >> p & 1) + 1 > caps.max4 for p in pts):
                continue
            for p in pts: used[p] += w
            fam.append(b)
            rec(i + 1, fam, used)
            fam.pop()
            for p in pts: used[p] -= w
    rec(0, [], [0] * n)
    return fams

def generator_families(n, caps):
    fams = phase1(n, list(range(5, n)), caps, verbose=False)
    out = set()
    for fam, aut in fams:
        if aut is None:
            ok, aut = is_canonical_marked(n, fam); assert ok
        ph = Phase2(n, fam, aut, caps, 0, collect_only=True)
        ph.run()
        for f in ph.collected:
            assert f not in out, f"duplicate {f}"
            out.add(f)
    return out

def V1():
    allok = True
    for n in (5, 6, 7):
        perms = all_perms(n)
        for mode in ('none', 'table'):
            caps = Caps(n, mode)
            t0 = time.time()
            lab = labelled_families(n, caps)
            canon = set(canon_bruteforce(n, f, perms) for f in lab)
            t1 = time.time()
            gen = generator_families(n, caps)
            gen_classes = set(canon_bruteforce(n, f, perms) for f in gen)
            t2 = time.time()
            ok = (canon == gen_classes) and len(gen_classes) == len(gen)   # same classes, no duplicates
            allok &= ok
            print(f"V1 n={n} mode={mode}: labelled families={len(lab)}, brute-force classes={len(canon)} "
                  f"[{t1-t0:.1f}s]; generator classes={len(gen)} [{t2-t1:.1f}s] => {'MATCH' if ok else 'MISMATCH'}",
                  flush=True)
            if not ok:
                print("   only brute:", sorted(canon - gen_classes)[:5]); print("   only gen:", sorted(gen_classes - canon)[:5])
    return allok

def V2():
    random.seed(1)
    allok = True
    for n in (7, 8):
        perms = all_perms(n)
        blocks = []
        for s in range(4, n):
            blocks += [list_to_mask(c) for c in itertools.combinations(range(n), s)]
        ncan = nfam = 0
        for trial in range(300):
            fam = []
            random.shuffle(blocks)
            for b in blocks[:random.randint(1, 8)]:
                if all(popcount(b & B) <= 2 for B in fam):
                    fam.append(b)
            fam = tuple(sorted(fam))
            cf = canon_bruteforce(n, fam, perms)
            for f in (fam, cf):
                ok, aut = is_canonical_marked(n, f)
                truth = (tuple(f) == cf)
                if ok != truth:
                    allok = False; print("V2 MISMATCH", n, f, cf, ok)
                if ok:
                    # automorphism group check by brute force
                    rows = np.sort(images(perms, list(f)), axis=1)
                    _, eq = rows_less_eq(rows, f)
                    bf = set(map(tuple, perms[eq].tolist())); mk = set(map(tuple, aut.tolist()))
                    if bf != mk or len(mk) != aut.shape[0]:
                        allok = False; print("V2 AUT MISMATCH", n, f, len(bf), len(mk), aut.shape)
                    ncan += 1
                nfam += 1
        print(f"V2 n={n}: {nfam} families tested, {ncan} canonical, automorphism groups verified => "
              f"{'OK' if allok else 'FAIL'}", flush=True)
    return allok

if __name__ == '__main__':
    ok1 = V1(); ok2 = V2()
    print("ALL OK" if ok1 and ok2 else "FAILURES")
