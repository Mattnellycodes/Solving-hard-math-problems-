"""Own hereditary Sylvester-Gallai check for candidate structures (c(9) audit).

Input: FAM lines from enum9 (or an explicit block list).  For each point p, the derived linear
space at p is {B - {p} : p in B}; after inversion about p these are exactly the lines with >= 3 points
of the real 8-point set P - {p} (the set is not collinear because P is not concyclic).
Sylvester-Gallai applied to a subset S of that 8-point set: if S is not collinear then some pair of S
spans a line containing no third point of S.  So any subset S (|S| >= 3) that is NOT inside one derived
line but in which EVERY pair lies in a derived line meeting S in >= 3 points is impossible.
The same test is applied to every maximum line set (lines at infinity) of the structure.
Also: lists the minimal violating subsets and identifies them (Fano plane = 7 points, 7 lines).
"""
import sys, itertools, re, ast

def parse(fname):
    out = []
    for line in open(fname):
        if line.startswith("FAM"):
            out.append([tuple(b) for b in ast.literal_eval(line.split("blocks=")[1].strip())])
    return out

def sg_closed_subsets(points, lines):
    """minimal subsets S of `points` violating SG w.r.t. `lines` (list of sets)."""
    pts = sorted(points)
    viol = []
    for r in range(3, len(pts) + 1):
        for S in itertools.combinations(pts, r):
            Sset = set(S)
            if any(Sset <= L for L in lines):
                continue
            if any(set(v) <= Sset for v in viol):
                continue
            ok = True
            for a, b in itertools.combinations(S, 2):
                if not any(a in L and b in L and len(L & Sset) >= 3 for L in lines):
                    ok = False; break
            if ok:
                viol.append(S)
    return viol

def max_line_sets(blocks, n, cap_lines):
    """all maximum sets of lines: lines = blocks or uncovered triples, pairwise sharing <= 1 point,
    pair coverage <= cap_lines."""
    cands = [set(b) for b in blocks]
    for t in itertools.combinations(range(n), 3):
        if not any(set(t) <= set(b) for b in blocks):
            cands.append(set(t))
    pairs = [len(c) * (len(c) - 1) // 2 for c in cands]
    best = [0]; sets = []
    def rec(start, chosen, used):
        if len(chosen) > best[0]:
            best[0] = len(chosen); sets.clear()
        if len(chosen) == best[0]:
            sets.append([sorted(cands[i]) for i in chosen])
        for j in range(start, len(cands)):
            if used + pairs[j] > cap_lines: continue
            if all(len(cands[j] & cands[i]) <= 1 for i in chosen):
                rec(j + 1, chosen + [j], used + pairs[j])
    rec(0, [], 0)
    return best[0], sets

def check(blocks, n=9, cap_lines=36):
    blocks = [set(b) for b in blocks]
    print("structure: sizes", sorted((len(b) for b in blocks), reverse=True), "blocks", [sorted(b) for b in blocks])
    killed = {}
    for p in range(n):
        derived = [b - {p} for b in blocks if p in b]
        v = sg_closed_subsets(set(range(n)) - {p}, derived)
        deg = len(derived)
        if v:
            killed[p] = v
            S = v[0]
            restricted = sorted(sorted(L & set(S)) for L in derived if len(L & set(S)) >= 3)
            print(f"  point {p} (degree {deg}): SG-closed subset {list(S)}; derived lines restricted: {restricted}"
                  f" -> {len(restricted)} lines on {len(S)} points"
                  + ("  [Fano plane]" if len(S) == 7 and len(restricted) == 7 and all(len(L) == 3 for L in restricted) else ""))
        else:
            print(f"  point {p} (degree {deg}): no SG-closed subset in the derived structure")
    ell, sets = max_line_sets([tuple(b) for b in blocks], n, cap_lines)
    bad = 0
    for L in sets:
        if sg_closed_subsets(set(range(n)), [set(l) for l in L]):
            bad += 1
    print(f"  maximum line sets: ellmax={ell}, {len(sets)} sets, {bad} of them violate hereditary SG at infinity")
    print(f"  => structure killed by hereditary SG at points {sorted(killed)}; killed at ALL points of degree>=7: "
          f"{all(p in killed for p in range(n) if sum(1 for b in blocks if p in b) >= 7)}")
    return bool(killed)

if __name__ == "__main__":
    seen = set()
    for fname in sys.argv[1:]:
        for blocks in parse(fname):
            key = tuple(sorted(tuple(sorted(b)) for b in blocks))
            if key in seen: continue
            seen.add(key)
            check(blocks)
