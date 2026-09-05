"""My own hereditary Sylvester-Gallai checker.  For a structure (n, blocks, lines):
 derived linear space at p = {B - {p} : p in B, |B| >= 4} on P - {p}  (plus, at infinity, the lines).
 A subset S (|S|>=3) is 'SG-closed' if S is not inside one derived line and every pair of S lies on a
 derived line containing >= 3 points of S.  By Sylvester-Gallai (applied to the inverted image of S)
 no SG-closed subset can exist in a real realisation."""
import itertools, json, sys
def sg_closed_subsets(points, lines, minimal_only=True):
    pts = sorted(points); found = []
    for r in range(3, len(pts) + 1):
        for S in itertools.combinations(pts, r):
            Ss = set(S)
            if any(Ss <= L for L in lines): continue
            if minimal_only and any(set(v) <= Ss for v in found): continue
            covered = set()
            for L in lines:
                I = L & Ss
                if len(I) >= 3:
                    covered |= set(itertools.combinations(sorted(I), 2))
            if len(covered) == r * (r - 1) // 2:
                found.append(S)
    return found
def check(n, blocks, lines):
    blocks = [frozenset(b) for b in blocks]; lines = [frozenset(l) for l in lines]
    out = {}
    for p in range(n):
        der = [b - {p} for b in blocks if p in b]
        v = sg_closed_subsets(set(range(n)) - {p}, der)
        if v: out[p] = v
    vinf = sg_closed_subsets(set(range(n)), lines)
    if vinf: out['inf'] = vinf
    return out
if __name__ == '__main__':
    recs = json.load(open(sys.argv[1]))
    for i, r in enumerate(recs):
        n = r.get('n', 9)
        v = check(n, r['blocks'], r['lines'])
        print(f'structure {i}: D+l={r["D_plus_l"]} circles={r["circles"]} sizes={sorted(len(b) for b in r["blocks"])} lines={len(r["lines"])} -> SG-closed subsets: {v if v else "NONE"}')
