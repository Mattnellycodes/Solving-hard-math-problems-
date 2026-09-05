"""Post-process enum9 output: isomorphism classes (invariant + exact canonical form), hereditary
Sylvester-Gallai check at every point (derived structure) written from scratch, and report."""
import sys, itertools, collections
from math import comb

def parse(fn):
    fams = []
    for line in open(fn):
        if not line.startswith("FAM"): continue
        parts = line.split()
        D = int(parts[1][2:]); ell = int(parts[2][4:]); cnt = int(parts[3][6:])
        blocks = [frozenset(int(ch) for ch in tok) for tok in parts[6:]]
        fams.append((D, ell, cnt, blocks))
    return fams

def invariant(blocks, n=9):
    deg = collections.Counter()
    for b in blocks:
        for p in b: deg[(p, len(b))] += 1
    prof = tuple(sorted(tuple(sorted((k[1], v) for k, v in deg.items() if k[0] == p)) for p in range(n)))
    pairs = collections.Counter()
    for b in blocks:
        for x, y in itertools.combinations(sorted(b), 2): pairs[(x, y)] += 1
    return (tuple(sorted(len(b) for b in blocks)), prof, tuple(sorted(pairs.values())))

def canon(blocks, n=9):
    best = None
    for perm in itertools.permutations(range(n)):
        img = tuple(sorted(tuple(sorted(perm[i] for i in b)) for b in blocks))
        if best is None or img < best: best = img
    return best

def sg_closed_subsets(points, lines):
    """Return a minimal subset S (|S|>=3, not inside one line) of `points` such that every pair of S
    lies on a line containing >= 3 points of S -- a Sylvester-Gallai violation -- or None."""
    pts = sorted(points)
    for r in range(3, len(pts) + 1):
        for S in itertools.combinations(pts, r):
            Sset = set(S)
            if any(Sset <= l for l in lines): continue
            covered = set()
            for l in lines:
                ls = l & Sset
                if len(ls) >= 3:
                    covered.update(itertools.combinations(sorted(ls), 2))
            if len(covered) == comb(r, 2):
                return S
    return None

def check_family(blocks, n=9):
    """hereditary SG at each point; returns dict point -> violating subset."""
    out = {}
    for p in range(n):
        derived = [b - {p} for b in blocks if p in b and len(b) >= 4]
        S = sg_closed_subsets(set(range(n)) - {p}, derived)
        if S: out[p] = S
    return out

if __name__ == "__main__":
    fams = parse(sys.argv[1])
    print("families read:", len(fams))
    byinv = collections.defaultdict(list)
    for f in fams: byinv[invariant(f[3])].append(f)
    print("distinct invariants:", len(byinv))
    classes = {}
    for inv, lst in byinv.items():
        # exact canonical forms within an invariant bucket
        cf = collections.defaultdict(list)
        for f in lst: cf[canon(f[3])].append(f)
        for c, l2 in cf.items(): classes[c] = l2
    print("isomorphism classes:", len(classes))
    for i, (c, l2) in enumerate(classes.items()):
        D, ell, cnt, blocks = l2[0]
        deg = [sum(1 for b in blocks if p in b) for p in range(9)]
        viol = check_family(blocks)
        print(f"class {i}: labelled copies={len(l2)} sizes={sorted((len(b) for b in blocks), reverse=True)} "
              f"D={D} ell_max={ell} count={cnt} degrees={deg}")
        print("   blocks:", [sorted(b) for b in blocks])
        print("   hereditary-SG violations at points:", {p: S for p, S in viol.items()})
        print("   -> killed" if viol else "   -> SURVIVES combinatorial filters")
