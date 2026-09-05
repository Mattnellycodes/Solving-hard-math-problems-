"""Skeptic c7-audit: brute-force enumeration of ALL labelled abstract structures (F, L) on 7 points
with C(7,3) - D - ell <= TARGET (default 10), written from scratch (no orderly generation, no
canonical-form pruning: every labelled family is visited, then classified up to isomorphism).

Conditions used:
  (C1) F = family of subsets of size 4..6, pairwise sharing <= 2 points;
       L = subset of F u {triples not inside any member of F}, pairwise sharing <= 1 point.
Cap regimes reported separately:
  none   : (C1) only
  weak   : Sylvester-Gallai only, o(m) >= 1 for all m:  derived cap 15-1 = 14 pairs at every point,
           line cap 21-1 = 20 pairs
  strong : o(6) = 3, o(7) = 3 (Kelly-Moser / classical): derived cap 12, line cap 18
Derived cap at p:  sum_{B in F, p in B} C(|B|-1, 2) <= cap.   Line cap: sum_{l in L} C(|l|,2) <= cap.
"""
import itertools, sys, json
from math import comb

n = 7
TARGET = int(sys.argv[1]) if len(sys.argv) > 1 else 10
N3 = comb(n, 3)
pts = list(range(n))
S = {k: [frozenset(c) for c in itertools.combinations(pts, k)] for k in (3, 4, 5, 6)}
cands = S[4] + S[5] + S[6]
deficit = {B: comb(len(B), 3) - 1 for B in cands}

families = []          # all labelled families (as tuples of frozensets)
def dfs(start, F):
    families.append(tuple(F))
    for j in range(start, len(cands)):
        B = cands[j]
        if all(len(B & A) <= 2 for A in F):
            F.append(B); dfs(j + 1, F); F.pop()
dfs(0, [])
print(f"labelled families F satisfying (C1): {len(families)}")

def line_cliques(F, need, line_cap):
    """all sets L of pairwise <=1-intersecting members of F u uncovered triples with |L| >= need
    and pair coverage <= line_cap; returns list of such L (all, not only maximal)."""
    covered = lambda t: any(t <= B for B in F)
    lc = list(F) + [t for t in S[3] if not covered(t)]
    out = []
    def rec(i, L, pairs):
        if len(L) >= need:
            out.append(tuple(L))
        for j in range(i, len(lc)):
            l = lc[j]
            if pairs + comb(len(l), 2) > line_cap: continue
            if all(len(l & m) <= 1 for m in L):
                L.append(l); rec(j + 1, L, pairs + comb(len(l), 2)); L.pop()
    rec(0, [], 0)
    return out

perms = list(itertools.permutations(pts))
def canon(F, L=()):
    best = None
    for s in perms:
        img = (tuple(sorted(tuple(sorted(s[i] for i in B)) for B in F)),
               tuple(sorted(tuple(sorted(s[i] for i in l)) for l in L)))
        if best is None or img < best: best = img
    return best

regimes = {"none": (10**9, 10**9), "weak": (14, 20), "strong": (12, 18)}
results = {r: {} for r in regimes}
n_with_D = 0
for F in families:
    D = sum(deficit[B] for B in F)
    need = N3 - TARGET - D
    if need > 7:   # at most 7 lines ever (each covers >= 3 of the 21 pairs)
        continue
    n_with_D += 1
    for r, (dcap, lcap) in regimes.items():
        if any(sum(comb(len(B) - 1, 2) for B in F if p in B) > dcap for p in pts):
            continue
        Ls = line_cliques(F, max(need, 0), lcap)
        for L in Ls:
            key = canon(F, L)
            results[r].setdefault(key, {"count": N3 - D - len(L), "D": D, "ell": len(L),
                                        "sizes": sorted((len(B) for B in F), reverse=True),
                                        "labelled": 0})
            results[r][key]["labelled"] += 1
print(f"families with D large enough to possibly reach target {TARGET}: {n_with_D}")
for r in regimes:
    print(f"\n=== regime {r}: {len(results[r])} isomorphism classes of (F, L) with circles <= {TARGET}")
    Fclasses = {}
    for key, rec in sorted(results[r].items(), key=lambda kv: (kv[1]['count'], kv[0])):
        Fclasses.setdefault(key[0], []).append(rec)
        print(f"  count={rec['count']} D={rec['D']} ell={rec['ell']} sizes={rec['sizes']} "
              f"blocks={[list(b) for b in key[0]]} lines={[list(l) for l in key[1]]} (#labelled {rec['labelled']})")
    print(f"  distinct block families F (up to iso): {len(Fclasses)}")
    for Fk in Fclasses:
        print("   F =", [list(b) for b in Fk], " ell values:", sorted(set(x['ell'] for x in Fclasses[Fk])))
json.dump({r: [{"F": [list(b) for b in k[0]], "L": [list(l) for l in k[1]], **v} for k, v in results[r].items()]
           for r in regimes}, open(f"enum_c7_t{TARGET}.json", "w"), indent=1)
