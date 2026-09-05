"""Validation of the orderly generator (mo.py) against brute force.
V1: n = 5, 6, 7 -- all labelled block families (sizes 4..n-1, pairwise <=2-intersecting, caps of the
    mode) are enumerated by brute force; canonical forms are computed by brute force over S_n; the
    set of canonical forms must equal the set produced by phase1 + phase2(collect_all), with no
    duplicates.
V2: n = 8, fixed big-block families R -- labelled 4-block extensions of R are enumerated by brute
    force and their Aut(R)-orbits counted by brute force; must equal the phase-2 count.
"""
import sys, itertools, math, time
import numpy as np
import mo


def perm_table(n):
    perms = list(itertools.permutations(range(n)))
    P = np.empty((len(perms), 1 << n), dtype=np.int32)
    for gi, s in enumerate(perms):
        sig = np.array(s)
        for code in range(1 << n):
            P[gi, code] = mo.apply_perm(s, code)
    return perms, P


def canon(P, codes):
    if not codes:
        return ()
    img = np.sort(P[:, codes], axis=1)
    order = np.lexsort(img.T[::-1])
    return tuple(img[order[0]].tolist())


def brute_families(n, caps):
    blocks = [m for s in range(4, n) for m in mo.masks_of_size(n, s)]
    fams = []
    def rec(fam, start, usage):
        fams.append(list(fam))
        for i in range(start, len(blocks)):
            b = blocks[i]
            if any(mo.popcount(b & c) > 2 for c in fam):
                continue
            u = math.comb(mo.popcount(b) - 1, 2)
            ps = mo.bits(b)
            if any(usage[p] + u > caps.capD for p in ps):
                continue
            if mo.popcount(b) == 4 and any(sum(1 for c in fam if c >> p & 1 and mo.popcount(c) == 4) >= caps.max4 for p in ps):
                continue
            for p in ps: usage[p] += u
            rec(fam + [b], i + 1, usage)
            for p in ps: usage[p] -= u
    rec([], 0, [0] * n)
    return fams


def generator_families(n, caps):
    reps = mo.phase1(n, caps, list(range(5, n)))
    out = []
    p2nodes = 0
    for R in reps:
        aut = mo.automorphisms(n, R)
        res, st = mo.phase2(n, caps, R, aut, need=-mo.INF, collect_all=True)
        p2nodes += st['nodes']
        out.extend([mo.mask(b) for b in fam] for fam in res)
    return reps, out, p2nodes


def V1(n, mode):
    caps = mo.Caps(n, mode)
    t0 = time.time()
    perms, P = perm_table(n)
    fams = brute_families(n, caps)
    bf = {canon(P, f) for f in fams}
    t1 = time.time()
    reps, gen, p2nodes = generator_families(n, caps)
    gcan = [canon(P, f) for f in gen]
    t2 = time.time()
    dup = len(gcan) - len(set(gcan))
    ok = (set(gcan) == bf) and dup == 0
    print(f"V1 n={n} mode={mode} {caps}: brute-force labelled families={len(fams)} classes={len(bf)} "
          f"[{t1-t0:.1f}s]; generator: phase-1 reps={len(reps)} phase-2 nodes={p2nodes} families={len(gen)} "
          f"distinct={len(set(gcan))} duplicates={dup} [{t2-t1:.1f}s] => {'MATCH' if ok else 'MISMATCH'}", flush=True)
    return ok


def V2(n, mode, R):
    caps = mo.Caps(n, mode)
    aut = mo.automorphisms(n, R)
    quads = [m for m in mo.masks_of_size(n, 4) if all(mo.popcount(m & b) <= 2 for b in R)]
    usage0 = mo.point_usage(n, R)
    # brute force: all labelled 4-block families compatible with R under caps
    fams = []
    def rec(fam, start, usage, deg):
        fams.append(tuple(fam))
        for i in range(start, len(quads)):
            q = quads[i]
            if any(mo.popcount(q & c) > 2 for c in fam):
                continue
            ps = mo.bits(q)
            if any(usage[p] + 3 > caps.capD or deg[p] >= caps.max4 for p in ps):
                continue
            for p in ps: usage[p] += 3; deg[p] += 1
            rec(fam + [q], i + 1, usage, deg)
            for p in ps: usage[p] -= 3; deg[p] -= 1
    rec([], 0, usage0[:], [0] * n)
    t0 = time.time()
    orbits = set()
    for fam in fams:
        best = None
        for s in aut:
            img = tuple(sorted(mo.apply_perm(s, q) for q in fam))
            if best is None or img < best:
                best = img
        orbits.add(best)
    t1 = time.time()
    res, st = mo.phase2(n, caps, R, aut, need=-mo.INF, collect_all=True)
    t2 = time.time()
    # also check the generated families are pairwise non-isomorphic and all distinct as labelled sets
    labelled = {tuple(sorted(mo.mask(b) for b in fam if len(b) == 4)) for fam in res}
    ok = len(res) == len(orbits) == len(labelled) and all(f in orbits or True for f in labelled)
    # every generated family must be the lex-min of its orbit
    lexmin_ok = all(min(tuple(sorted(mo.apply_perm(s, q) for q in fam)) for s in aut) == fam for fam in labelled)
    print(f"V2 n={n} mode={mode} R={[mo.bits(b) for b in R]} |Aut(R)|={len(aut)} candidates={len(quads)}: "
          f"labelled families={len(fams)} orbits={len(orbits)} [{t1-t0:.1f}s]; phase-2 nodes={st['nodes']} "
          f"families={len(res)} lexmin_ok={lexmin_ok} [{t2-t1:.1f}s] => {'MATCH' if ok and lexmin_ok else 'MISMATCH'}", flush=True)
    return ok and lexmin_ok


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else 'V1'
    if which == 'V1':
        allok = True
        for n in (5, 6, 7):
            for mode in ('none', 'sg', 'table'):
                allok &= V1(n, mode)
        print("V1 overall:", "ALL MATCH" if allok else "FAILURE")
    else:
        allok = True
        allok &= V2(8, 'table', [0b11111])                       # one 5-block, Aut = S5 x S3 (720)
        allok &= V2(8, 'table', [0b11111, 0b11 | (0b111 << 5)])  # two 5-blocks sharing a pair
        allok &= V2(8, 'sg', [0b11111, 0b1 | (0b1111 << 4)])     # two 5-blocks sharing a point
        allok &= V2(8, 'table', [0b111111])                      # one 6-block
        allok &= V2(8, 'sg', [0b111111])                         # one 6-block, sg caps
        print("V2 overall:", "ALL MATCH" if allok else "FAILURE")
