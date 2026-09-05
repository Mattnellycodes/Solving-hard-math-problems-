"""Necessary realisability conditions for an abstract structure (F, L) on P = {0..9}
(F = rich blocks of size >= 4, L = lines of size >= 3), written independently of ../theory and
../verify_independent.

1. hereditary Sylvester–Gallai.  Inverting about p in P turns the blocks through p into the >= 3-point
   lines of a real 9-point set Q_p = P \ {p} (the 'derived structure at p'); the lines L are the
   >= 3-point lines of the real 10-point set P itself.  For any real point set X and any S ⊆ X with
   |S| >= 3 not collinear, S has an ordinary line, i.e. a pair of S on no >= 3-point line of X that
   contains a third point of S.  A violating S is 'SG-closed'.  We test every S (2^9 resp. 2^10 subsets).
2. no derived (8_3): in a real point set, no 8 points carry 8 three-point lines
   (Möbius–Kantor, proved unrealisable in ../verify_independent/mk_unrealisable2.py); we test every
   8-subset S of Q_p (resp. of P) and count lines meeting S in exactly 3 points.  (7 lines on 7 points is
   Fano and is caught by 1., but we count it too.)
3. Miquel closure on the 11-point Möbius set P ∪ {∞}: blocks are F together with S ∪ {∞} for S in L.
   Miquel's theorem (real inversive plane): if 8 distinct points are labelled by the vertices of a cube
   and 5 of the 6 faces are concyclic quadruples, so is the 6th.  In an abstract structure a quadruple is
   concyclic iff it is contained in a block; a labelled cube with exactly 5 concyclic faces is a violation.
"""
import itertools
from collections import Counter

INF = 10


def derived(F, p):
    return [frozenset(B - {p}) for B in F if p in B]


def sg_closed(points, lines):
    """Return the list of SG-closed subsets (minimal ones only) of the linear space (points, lines)."""
    pts = sorted(points)
    lines = [frozenset(l) for l in lines]
    pair_line = {}
    for l in lines:
        for a, b in itertools.combinations(sorted(l), 2):
            pair_line[(a, b)] = l
    viol = []
    for r in range(3, len(pts) + 1):
        for S in itertools.combinations(pts, r):
            Sset = frozenset(S)
            if any(Sset <= l for l in lines):
                continue
            ok = False
            for a, b in itertools.combinations(S, 2):
                l = pair_line.get((a, b))
                if l is None or len(l & Sset) == 2:
                    ok = True; break
            if not ok and not any(v <= Sset for v in viol):
                viol.append(Sset)
    return viol


def hereditary_sg(F, L, P=range(10)):
    """dict: point -> violating subsets (derived structures), key 'inf' for the line structure."""
    out = {}
    for p in P:
        v = sg_closed(set(P) - {p}, derived(F, p))
        if v:
            out[p] = v
    v = sg_closed(set(P), L)
    if v:
        out['inf'] = v
    return out


def three_line_excess(F, L, P=range(10)):
    """Violations of t3(7) <= 6 (Fano) and t3(8) <= 7 ((8_3)) in derived structures and lines."""
    out = []
    def check(points, lines, tag):
        for r, cap in ((7, 6), (8, 7)):
            for S in itertools.combinations(sorted(points), r):
                Sset = frozenset(S)
                c = sum(1 for l in lines if len(l & Sset) == 3)
                if c > cap:
                    out.append((tag, r, S, c))
    for p in P:
        check(set(P) - {p}, derived(F, p), ('derived', p))
    check(set(P), [frozenset(l) for l in L], ('lines',))
    return out


def three_line_counts(F, L, P=range(10)):
    d = {p: sum(1 for B in F if p in B and len(B) == 4) for p in P}
    return d, sum(1 for S in L if len(S) == 3)


_CUBES = None
def cube_labellings():
    """All 840 distinct face-systems of a cube on vertex labels 0..7."""
    global _CUBES
    if _CUBES is None:
        verts = list(itertools.product((0, 1), repeat=3))
        faces = [[v for v in verts if v[i] == b] for i in range(3) for b in (0, 1)]
        seen = set()
        for perm in itertools.permutations(range(8)):
            fs = frozenset(frozenset(perm[verts.index(v)] for v in f) for f in faces)
            seen.add(fs)
        _CUBES = [sorted(tuple(sorted(f)) for f in fs) for fs in seen]
    return _CUBES


def miquel_violations(F, L, P=range(10)):
    blocks = [frozenset(B) for B in F] + [frozenset(S) | {INF} for S in L]
    pts = list(P) + [INF]
    def coblock(Q):
        return any(Q <= B for B in blocks)
    quads = {}
    for Q in itertools.combinations(pts, 4):
        quads[frozenset(Q)] = coblock(frozenset(Q))
    viol = []
    cubes = cube_labellings()
    for S in itertools.combinations(pts, 8):
        for fs in cubes:
            faces = [frozenset(S[i] for i in f) for f in fs]
            c = sum(quads[f] for f in faces)
            if c == 5:
                viol.append((S, [tuple(sorted(f)) for f in faces if not quads[f]][0]))
    return viol


def full_report(F, L, verbose=True):
    F = [frozenset(B) for B in F]; L = [frozenset(S) for S in L]
    hs = hereditary_sg(F, L)
    tl = three_line_excess(F, L)
    mq = miquel_violations(F, L)
    d4, l3 = three_line_counts(F, L)
    rep = {'hereditary_sg': {str(k): [sorted(s) for s in v] for k, v in hs.items()},
           'three_line_excess': tl, 'miquel': mq[:5], 'n_miquel': len(mq), 'd4': d4, 'l3': l3}
    if verbose:
        print(f"   hereditarySG violations at: {sorted(hs.keys(), key=str)}  |  Fano/(8_3) excess: {len(tl)}  |  "
              f"Miquel violations: {len(mq)}  |  d4={[d4[p] for p in range(10)]} three-lines={l3}")
    return rep


if __name__ == '__main__':
    # self-tests
    FANO = [{0,1,2},{0,3,4},{0,5,6},{1,3,5},{1,4,6},{2,3,6},{2,4,5}]
    assert sg_closed(range(7), FANO) == [frozenset(range(7))]
    MK = [{i % 8, (i+1) % 8, (i+3) % 8} for i in range(8)]
    assert sg_closed(range(8), MK) == []          # MK has ordinary lines, not caught by SG
    assert len(cube_labellings()) == 840
    # cube structure on 8 points (SQS(8) minus a parallel class) is Miquel-closed: 12 blocks + 8 uncovered
    verts = list(itertools.product((0, 1), repeat=3))
    idx = {v: i for i, v in enumerate(verts)}
    faces = [frozenset(idx[v] for v in verts if v[i] == b) for i in range(3) for b in (0, 1)]
    diag = []
    for i, j in itertools.combinations(range(3), 2):
        for s in (0, 1):
            diag.append(frozenset(idx[v] for v in verts if (v[i] + v[j]) % 2 == s))
    cube = faces + diag
    # embed in 10 points: F = cube blocks on points 0..7, no lines
    print("cube structure Miquel violations:", len(miquel_violations(cube, [])))
    # remove one face: 5 faces present -> violation expected
    print("cube minus a face Miquel violations:", len(miquel_violations(cube[1:], [])))
    print("self-tests passed")


def bundle_violations(F, L, P=range(10)):
    """Bundle theorem on P ∪ {∞}: four pairwise disjoint pairs with exactly 5 of the 6 unions in blocks."""
    blocks = [frozenset(B) for B in F] + [frozenset(S) | {INF} for S in L]
    pts = list(P) + [INF]
    quads = {}
    def cob(Q):
        if Q not in quads:
            quads[Q] = any(Q <= B for B in blocks)
        return quads[Q]
    pairs = [frozenset(t) for t in itertools.combinations(pts, 2)]
    viol = []
    def rec(chosen, start, used):
        if len(chosen) == 4:
            c = sum(cob(Pi | Pj) for Pi, Pj in itertools.combinations(chosen, 2))
            if c == 5:
                viol.append([tuple(sorted(p)) for p in chosen])
            return
        for i in range(start, len(pairs)):
            if pairs[i] & used:
                continue
            rec(chosen + [pairs[i]], i + 1, used | pairs[i])
    rec([], 0, frozenset())
    return viol


def fregier_violations(F, L, P=range(10)):
    """Two-block Frégier lemma (sphere model), applied to the 11-point Möbius set P ∪ {∞}.
    Let B0, B1 be distinct blocks (circles, i.e. plane sections Π0, Π1 of the sphere; B1 may be an
    uncovered triple), m = Π0 ∩ Π1.  For pairs A ⊂ B0 \ B1 and Q ⊂ B1 \ B0 whose union lies in a rich
    block K: the planes Π0, Π1, Π_K meet pairwise in line(A), line(Q), m, which are concurrent, so
    line(A) ∩ m = line(Q) ∩ m.  (line(A) ≠ m since A ⊄ B1.)  In the bipartite graph 'inside pairs of
    B0 \ B1' vs 'pairs of B1 \ B0' with an edge when the union is in a rich block, all chord-lines of a
    connected component pass through one point X ∈ m, and X is not on the sphere (X ∈ line(Q) whose
    sphere points are Q ⊄ B0).  Hence inside pairs of a component are pairwise disjoint (two chords of
    circle(B0) through a common point of B0 would put X on the sphere) and outside pairs are pairwise
    disjoint (two secants of the sphere through a common sphere point and X would coincide)."""
    rich = [frozenset(B) for B in F] + [frozenset(S) | {INF} for S in L]
    pts = list(P) + [INF]
    # all blocks of the Möbius set: rich blocks and uncovered triples
    covered = set()
    for K in rich:
        for T in itertools.combinations(sorted(K), 3):
            covered.add(frozenset(T))
    blocks = rich + [frozenset(T) for T in itertools.combinations(pts, 3) if frozenset(T) not in covered]
    richset = rich
    def coblocked(Q4):
        return any(Q4 <= K for K in richset)
    viol = []
    for B0, B1 in itertools.permutations(blocks, 2):
        D0 = sorted(B0 - B1); D1 = sorted(B1 - B0)
        if len(D0) < 2 or len(D1) < 2:
            continue
        inside = [frozenset(t) for t in itertools.combinations(D0, 2)]
        outside = [frozenset(t) for t in itertools.combinations(D1, 2)]
        adj = {}
        for A in inside:
            for Q in outside:
                if coblocked(A | Q):
                    adj.setdefault(A, set()).add(Q); adj.setdefault(Q, set()).add(A)
        seen = set()
        for v in inside:
            if v in seen or v not in adj:
                continue
            comp = set(); stack = [v]
            while stack:
                u = stack.pop()
                if u in comp:
                    continue
                comp.add(u); stack.extend(adj[u])
            seen |= comp
            ins = [u for u in comp if u in set(inside)]; outs = [u for u in comp if u in set(outside)]
            if any(a & b for a, b in itertools.combinations(ins, 2)) or any(a & b for a, b in itertools.combinations(outs, 2)):
                viol.append((tuple(sorted(B0)), tuple(sorted(B1)), [tuple(sorted(u)) for u in ins], [tuple(sorted(u)) for u in outs]))
    return viol


def derived11(F, L, p, P=range(10)):
    """Derived structure at p in the 11-point Möbius set P ∪ {∞}: a real 10-point set (P \\ {p}) ∪ {∞'}
    whose >= 3-point lines are B - {p} for rich blocks B ∋ p and (S - {p}) ∪ {∞'} for lines S ∋ p."""
    pts = set(P) - {p} | {INF}
    lines = [frozenset(B) - {p} for B in F if p in B] + [(frozenset(S) - {p}) | {INF} for S in L if p in S]
    return pts, lines


def eleven_point_report(F, L, P=range(10)):
    """For every p: pairs covered by the derived 10-point structure (must be <= 45 - o(10); SG: <= 44,
    Kelly–Moser/Csima–Sawyer: <= 40), number of 3-point lines (<= t3(10) = 12, proved in t3_10_proof.py),
    SG-closed subsets, Fano/(8_3) excess.  Returns dict of violations under SG-only assumptions and the
    max pair count (to flag structures excluded only by o(10) >= 5)."""
    F = [frozenset(B) for B in F]; L = [frozenset(S) for S in L]
    out = {'sg_closed': {}, 'three_line_excess': [], 'orchard10': [], 'pairs': {}}
    for p in P:
        pts, lines = derived11(F, L, p)
        pairs = sum(len(l) * (len(l) - 1) // 2 for l in lines)
        out['pairs'][p] = pairs
        n3 = sum(1 for l in lines if len(l) == 3)
        if n3 > 12:
            out['orchard10'].append((p, n3))
        v = sg_closed(pts, lines)
        if v:
            out['sg_closed'][p] = [sorted(x, key=str) for x in v]
        for r, cap in ((7, 6), (8, 7)):
            for S_ in itertools.combinations(sorted(pts, key=str), r):
                Sset = frozenset(S_)
                c = sum(1 for l in lines if len(l & Sset) == 3)
                if c > cap:
                    out['three_line_excess'].append((p, r, c))
    out['max_pairs'] = max(out['pairs'].values())
    out['dead_sg_only'] = bool(out['sg_closed']) or bool(out['three_line_excess']) or bool(out['orchard10']) or out['max_pairs'] > 44
    out['dead_o10_5'] = out['dead_sg_only'] or out['max_pairs'] > 40
    return out
