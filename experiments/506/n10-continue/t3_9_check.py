"""Own verification of t3(9) <= 10 given the (8_3) lemma: every family of 11 triples on 9 points pairwise
sharing <= 1 point contains 8 points carrying 8 of the triples (a Möbius–Kantor (8_3)), which is not
realisable by 8 distinct real points; hence 9 real points carry <= 10 three-point lines.
Counting: 11 triples cover 33 of the 36 pairs; every point has even leave-degree (8 - 2k), so the 3 uncovered
pairs form a triangle; up to relabelling the leave is {6,7,8}.  All decompositions of K9 minus that triangle
into 11 triangles are enumerated by DFS."""
import itertools
LEAVE = {frozenset((6, 7)), frozenset((6, 8)), frozenset((7, 8))}
EDGES = [frozenset(e) for e in itertools.combinations(range(9), 2) if frozenset(e) not in LEAVE]
fams = []
def dfs(uncov, chosen):
    if not uncov:
        fams.append(chosen[:]); return
    e = min(uncov, key=lambda f: tuple(sorted(f))); a, b = sorted(e)
    for c in range(9):
        if c in e: continue
        if frozenset((a, c)) in uncov and frozenset((b, c)) in uncov:
            chosen.append(frozenset((a, b, c)))
            dfs(uncov - {e, frozenset((a, c)), frozenset((b, c))}, chosen)
            chosen.pop()
dfs(frozenset(EDGES), [])
# also check that the parity argument is right: leave of any 11-packing is a triangle (brute force over all packings)
def all_packings():
    out = []
    def rec(start, chosen, used):
        if len(chosen) == 11:
            out.append(chosen[:]); return
        for i in range(start, len(TR)):
            T = TR[i]
            if not (PAIRS[i] & used):
                chosen.append(T); rec(i + 1, chosen, used | PAIRS[i]); chosen.pop()
    rec(0, [], frozenset())
    return out
TR = [frozenset(t) for t in itertools.combinations(range(9), 3)]
PAIRS = [frozenset(frozenset(p) for p in itertools.combinations(sorted(T), 2)) for T in TR]
allp = all_packings()
leaves_tri = all(len({v for e in (set(frozenset(p) for p in itertools.combinations(range(9), 2)) - frozenset().union(*[PAIRS[TR.index(T)] for T in fam])) for v in e}) == 3 for fam in allp)
has83 = all(any(sum(1 for T in fam if T <= set(sub)) == 8 for sub in itertools.combinations(range(9), 8)) for fam in fams)
print(f"11-packings with leave triangle {{6,7,8}}: {len(fams)}; all labelled 11-packings: {len(allp)} (= 84 triangles x {len(fams)}: {len(allp) == 84 * len(fams)}); every leave is a triangle: {leaves_tri}")
print(f"every 11-packing contains an (8_3): {has83}  => t3(9) <= 10 given the (8_3) lemma")
