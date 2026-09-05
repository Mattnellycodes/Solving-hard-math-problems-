"""lemmas-audit: post-process enum_own output: isomorphism classes, line sets, geometric filters.
usage: python3 post_own.py n capmode file.out [--lines]
Own implementation of:
  * isomorphism classification (networkx VF2 on the incidence bipartite graph, blocks coloured by size)
  * ell_max and all maximum line sets (pairwise <=1 sharing, pair budget)
  * hereditary Sylvester-Gallai test on every derived structure (blocks through p, minus p) and on the
    line sets: a subset S (|S|>=3) NOT contained in a single (derived) line such that every pair of S is
    covered by a (derived) line meeting S in >= 3 points is a violation (S inverted would be a
    non-collinear real set with no ordinary line).
  * degree-8 test (derived (8_3)), by Lemma 2.4.
"""
import sys, itertools
from math import comb
import networkx as nx

def olower(m, capmode):
    if m <= 2: return 0
    if capmode == 0: return 0
    if capmode == 1: return 1
    tab = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5}
    return tab[m] if m in tab else -(-6 * m // 13)

def decode(mask, n):
    return frozenset(i for i in range(n) if mask >> i & 1)

def graph_of(family, n):
    G = nx.Graph()
    for p in range(n):
        G.add_node(("p", p), kind="p")
    for i, b in enumerate(family):
        G.add_node(("b", i), kind=f"b{len(b)}")
        for p in b:
            G.add_edge(("b", i), ("p", p))
    return G

def invariant(family, n):
    degs = sorted(sum(1 for b in family if p in b) for p in range(n))
    sizes = sorted(len(b) for b in family)
    pairc = sorted(sum(1 for b in family if a in b and c in b) for a, c in itertools.combinations(range(n), 2))
    return (tuple(sizes), tuple(degs), tuple(pairc))

def sg_violations(points, lines):
    """minimal subsets S violating hereditary SG (returns at most a few)."""
    pts = sorted(points)
    viol = []
    for r in range(3, len(pts) + 1):
        for S in itertools.combinations(pts, r):
            Ss = frozenset(S)
            if any(Ss <= l for l in lines):
                continue
            covered = set()
            for l in lines:
                L = l & Ss
                if len(L) >= 3:
                    covered.update(itertools.combinations(sorted(L), 2))
            if len(covered) == comb(r, 2):
                if not any(frozenset(v) <= Ss for v in viol):
                    viol.append(S)
        if len(viol) >= 3:
            break
    return viol


def complete_quadrilateral_pattern(points, lines):
    """lines: 4 triples on 6 points, pairwise sharing exactly 1 point, each point on exactly 2 lines.
    Returns the three 'opposite pairs' (pairs of points on no common line) or None."""
    pts = sorted(points)
    if len(pts) != 6 or len(lines) != 4 or any(len(l) != 3 for l in lines):
        return None
    if any(len(a & b) != 1 for a, b in itertools.combinations(lines, 2)):
        return None
    if any(sum(1 for l in lines if p in l) != 2 for p in pts):
        return None
    opp = [frozenset(pr) for pr in itertools.combinations(pts, 2) if not any(set(pr) <= l for l in lines)]
    return opp if len(opp) == 3 else None

def angle_lemma_kills(fam, n, p, derived):
    """Angle-Lemma pattern at p: derived lines form a complete quadrilateral on Q = P\\{p} and the
    blocks avoiding p are exactly the three complements (in Q) of the opposite pairs."""
    Q = set(range(n)) - {p}
    opp = complete_quadrilateral_pattern(Q, derived)
    if opp is None:
        return False
    avoiding = [b for b in fam if p not in b]
    targets = {frozenset(Q - o) for o in opp}
    return set(avoiding) == targets and len(avoiding) == 3

def mk_at_infinity(L, n):
    """8 lines of size 3 with union of size 8, pairwise <= 1, and no further line of L meeting that
    8-set in >= 3 points: the rich lines of those 8 points form an (8_3) -> Lemma 2.4."""
    tri = [l for l in L if len(l) == 3]
    for comb8 in itertools.combinations(tri, 8):
        S = frozenset().union(*comb8)
        if len(S) != 8:
            continue
        others = [l for l in L if l not in comb8 and len(l & S) >= 3]
        if not others:
            return True
    return False

def line_sets(family, n, capl):
    cands = list(family) + [frozenset(t) for t in itertools.combinations(range(n), 3)
                             if not any(frozenset(t) <= b for b in family)]
    best = [0]; sets = []
    def rec(start, chosen, pairs):
        if len(chosen) > best[0]:
            best[0] = len(chosen); sets.clear()
        if len(chosen) == best[0]:
            sets.append(list(chosen))
        avail = [j for j in range(start, len(cands)) if all(len(cands[j] & c) <= 1 for c in chosen)
                 and pairs + comb(len(cands[j]), 2) <= capl]
        if len(chosen) + len(avail) < best[0]:
            return
        for j in avail:
            rec(j + 1, chosen + [cands[j]], pairs + comb(len(cands[j]), 2))
    rec(0, [], 0)
    return best[0], sets

def main():
    n = int(sys.argv[1]); capmode = int(sys.argv[2]); fn = sys.argv[3]
    capl = comb(n, 2) - olower(n, capmode)
    recs = []
    for line in open(fn):
        if line.startswith("REC"):
            head, masks = line.split(":")
            _, D, e, count = head.split()
            fam = [decode(int(m), n) for m in masks.split()]
            recs.append((int(D), int(e), int(count), fam))
        elif line.startswith("DONE"):
            print(line.strip())
    print(f"{len(recs)} labelled records")
    classes = []   # (inv, graph, representative, count)
    for D, e, count, fam in recs:
        inv = invariant(fam, n); G = graph_of(fam, n)
        for cl in classes:
            if cl[0] == inv and nx.is_isomorphic(cl[1], G, node_match=lambda a, b: a["kind"] == b["kind"]):
                cl[3] += 1; break
        else:
            classes.append([inv, G, (D, e, count, fam), 1])
    print(f"{len(classes)} isomorphism classes")
    for inv, G, (D, e, count, fam), mult in classes:
        sizes = sorted((len(b) for b in fam), reverse=True)
        degs = [sum(1 for b in fam if p in b) for p in range(n)]
        print(f"\nCLASS: sizes={sizes} D={D} ell_max={e} count={count} degrees={degs} labelled copies={mult}")
        print("  blocks:", [sorted(b) for b in fam])
        # own ell_max recomputation and line sets
        e2, sets = line_sets(fam, n, capl)
        print(f"  own ell_max={e2} (enumerator said {e}); #maximum line sets={len(sets)}; example={[sorted(l) for l in sets[0]] if sets else None}")
        # derived structures
        killed = {}
        for p in range(n):
            derived = [b - {p} for b in fam if p in b]
            if len(derived) >= 8 and all(len(d) == 3 for d in derived):
                killed[p] = "degree-8 point: derived structure is an (8_3) -> Lemma 2.4"
                continue
            v = sg_violations(set(range(n)) - {p}, derived)
            if v:
                killed[p] = f"hereditary SG violated on S={v[0]} (derived lines {[sorted(d) for d in derived]})"
        for p in range(n):
            derived = [b - {p} for b in fam if p in b]
            if p not in killed and angle_lemma_kills(fam, n, p, derived):
                killed[p] = "Angle Lemma: derived lines = complete quadrilateral, blocks avoiding p = its three diagonal quadruples"
        if killed:
            for p, why in killed.items():
                print(f"  point {p}: {why}")
        else:
            print("  no derived-structure violation (structure survives block-level filters)")
        if sets:
            # line-level filters on every maximum line set
            mk = sum(1 for L in sets if mk_at_infinity(L, n))
            cq = 0
            for L in sets:
                opp = complete_quadrilateral_pattern(set(range(n)), L) if n == 6 else None
                if opp is not None and {frozenset(set(range(n)) - o) for o in opp} == set(fam):
                    cq += 1
            print(f"  maximum line sets: {len(sets)}; killed by (8_3)-at-infinity: {mk}; killed by Angle Lemma at infinity (n=6 pattern): {cq}")
        # lines: check every maximum line set (and report if all violate)
        if "--lines" in sys.argv and sets:
            bad = 0
            for L in sets:
                if sg_violations(set(range(n)), L):
                    bad += 1
            print(f"  maximum line sets violating hereditary SG: {bad} of {len(sets)}")

if __name__ == "__main__":
    main()
