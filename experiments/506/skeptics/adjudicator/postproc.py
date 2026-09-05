"""Post-process FOUND lines of enum_big: isomorphism classes (networkx VF2 on the point-block
incidence graph), degrees, own hereditary-SG test (derived Fano: at a point p, a 7-subset S of
P\\{p} not inside one block through p, all of whose 21 pairs lie on derived lines with >= 3
points of S), own (8_3) test (a point in >= 8 four-blocks and no 5-block through it => derived
structure is an (8_3) => impossible), and own Angle-Lemma pattern test at points of degree 4
(4 derived 3-lines forming a complete quadrilateral whose 3 diagonal quadruples are blocks).
Usage: python3 postproc.py file [file ...]
"""
import sys, re
from itertools import combinations
from math import comb
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher


def parse(fn):
    fams = []
    for line in open(fn):
        if not line.startswith('FOUND'):
            continue
        D = int(re.search(r'D=(\d+)', line).group(1)); e = int(re.search(r'ell=(\d+)', line).group(1))
        cnt = int(re.search(r'count=(-?\d+)', line).group(1))
        blocks = [frozenset(int(x) for x in b.split(',')) for b in re.findall(r'\[([0-9,]+)\]', line)]
        fams.append((D, e, cnt, blocks))
    return fams


def graph(blocks, n):
    G = nx.Graph()
    for p in range(n):
        G.add_node(('p', p), kind='p')
    for i, B in enumerate(blocks):
        G.add_node(('b', i), kind='b%d' % len(B))
        for p in B:
            G.add_edge(('p', p), ('b', i))
    return G


def derived_fano(blocks, n):
    """return list of (p, S) where hereditary SG is violated by a derived Fano plane"""
    out = []
    for p in range(n):
        thr = [B - {p} for B in blocks if p in B]
        others = [q for q in range(n) if q != p]
        for S in combinations(others, 7):
            S = frozenset(S)
            if any(S <= L for L in thr):
                continue
            lines = [L & S for L in thr if len(L & S) >= 3]
            pairs = set()
            for L in lines:
                pairs |= {frozenset(pr) for pr in combinations(L, 2)}
            if len(pairs) == 21:
                out.append((p, sorted(S)))
                break
    return out


def mk_points(blocks, n):
    """points p with >= 8 four-blocks through p and no block of size >= 5 through p"""
    out = []
    for p in range(n):
        thr = [B for B in blocks if p in B]
        if len(thr) >= 8 and all(len(B) == 4 for B in thr):
            out.append((p, len(thr)))
    return out


def angle_pattern(blocks, n):
    """points p of degree exactly 4 (all 4-blocks) whose derived lines form a complete quadrilateral
    on P\\{p} (n=7: 6 other points) with the 3 blocks avoiding p equal to the 3 diagonal quadruples"""
    out = []
    if n != 7:
        return out
    for p in range(n):
        thr = [B - {p} for B in blocks if p in B]
        avoid = [B for B in blocks if p not in B]
        if len(thr) != 4 or any(len(L) != 3 for L in thr) or len(avoid) != 3:
            continue
        if not all(len(A & B) == 1 for A, B in combinations(thr, 2)):
            continue
        V = {(i, j): next(iter(thr[i] & thr[j])) for i, j in combinations(range(4), 2)}
        opp = [{V[(0, 1)], V[(2, 3)]}, {V[(0, 2)], V[(1, 3)]}, {V[(0, 3)], V[(1, 2)]}]
        diag = {frozenset(opp[0] | opp[1]), frozenset(opp[0] | opp[2]), frozenset(opp[1] | opp[2])}
        if diag == set(avoid):
            out.append(p)
    return out


for fn in sys.argv[1:]:
    fams = parse(fn)
    n = max((max(B) for _, _, _, bl in fams for B in bl), default=-1) + 1
    print(f"== {fn}: {len(fams)} labelled families, n={n}")
    reps = []  # (graph, info, count)
    for D, e, cnt, bl in fams:
        G = graph(bl, n)
        for r in reps:
            if r[2][0] == D and r[2][1] == e and GraphMatcher(r[0], G, node_match=lambda a, b: a['kind'] == b['kind']).is_isomorphic():
                r[3] += 1
                break
        else:
            reps.append([G, bl, (D, e, cnt), 1])
    print(f"   isomorphism classes: {len(reps)}")
    for G, bl, (D, e, cnt), mult in reps:
        degs = [sum(1 for B in bl if p in B) for p in range(n)]
        sizes = sorted(len(B) for B in bl)
        print(f"   CLASS x{mult}: sizes={sizes} D={D} ell_max={e} count={cnt} degrees={degs}")
        print(f"      blocks={[sorted(B) for B in bl]}")
        df = derived_fano(bl, n)
        if df:
            print(f"      hereditary SG violated (derived Fano) at points {[p for p, S in df]} e.g. p={df[0][0]} S={df[0][1]}")
        mk = mk_points(bl, n)
        if mk:
            print(f"      points in >= 8 four-blocks (derived (8_3)): {mk}")
        ap = angle_pattern(bl, n)
        if ap:
            print(f"      Angle-Lemma pattern (complete quadrilateral + 3 diagonal blocks) at points {ap}")
        if not df and not mk and not ap:
            print("      NOT killed by derived Fano / (8_3) / Angle pattern at block level")
