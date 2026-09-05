"""Audit toolkit for Erdős #506, n = 10 (written from scratch by the auditor).

A *structure* on P = {0..n-1} is (F, L): F = list of rich blocks (subsets of size >= 4, pairwise
sharing <= 2 points), L = list of lines (subsets of size >= 3; a line of size >= 4 must be a member
of F; a line of size 3 must be a triple not contained in any member of F; lines pairwise share <= 1
point).  D = sum_{B in F} (C(|B|,3) - 1), count = C(n,3) - D - |L|.

Real linear spaces forced by (F, L) on P (all are sets of DISTINCT real points, so Sylvester-Gallai,
the (8_3) lemma and the orchard caps apply):
  * derived 9-set at p:  points P - p, lines {B - p : B in F, p in B}
  * Q_p (10 points):     points (P - p) u {oo'}, lines {B - p (+oo' if B in L) : B in F, p in B}
                         u {(S - p) u {oo'} : S in L, |S| = 3, p in S}
  * (P, L) itself.
"""
import itertools
from math import comb
import networkx as nx
from networkx.algorithms import isomorphism

N = 10
P = tuple(range(N))


def pairs_of(s):
    return itertools.combinations(sorted(s), 2)


def D_of(F):
    return sum(comb(len(B), 3) - 1 for B in F)


def count_of(F, L, n=N):
    return comb(n, 3) - D_of(F) - len(L)


# ----------------------------------------------------------------------------------------------
# packing checks
def check_structure(F, L, n=N):
    """Raise AssertionError if (F, L) violates the packing axioms."""
    F = [frozenset(B) for B in F]
    L = [frozenset(S) for S in L]
    for B in F:
        assert len(B) >= 4
    for A, B in itertools.combinations(F, 2):
        assert len(A & B) <= 2, (sorted(A), sorted(B))
    Fset = set(F)
    for S in L:
        assert len(S) >= 3
        if len(S) >= 4:
            assert S in Fset, sorted(S)
        else:
            assert not any(S <= B for B in F), sorted(S)
    for S, T in itertools.combinations(L, 2):
        assert len(S & T) <= 1, (sorted(S), sorted(T))
    return True


# ----------------------------------------------------------------------------------------------
# canonical forms via coloured incidence graphs (exact isomorphism with VF2, bucketed by WL hash)
def incidence_graph(F, L, n=N):
    G = nx.Graph()
    for p in range(n):
        G.add_node(('p', p), c=0)
    Lset = set(frozenset(S) for S in L)
    for i, B in enumerate(F):
        col = 2 if frozenset(B) in Lset else 1
        G.add_node(('b', i), c=col)
        for p in B:
            G.add_edge(('b', i), ('p', p))
    for j, S in enumerate(L):
        if len(S) == 3:
            G.add_node(('l', j), c=3)
            for p in S:
                G.add_edge(('l', j), ('p', p))
    return G


def wl_hash(G):
    return nx.weisfeiler_lehman_graph_hash(G, node_attr='c', iterations=4)


def isomorphic(G1, G2):
    nm = isomorphism.categorical_node_match('c', -1)
    return isomorphism.GraphMatcher(G1, G2, node_match=nm).is_isomorphic()


class ClassCollector:
    """Collect structures up to isomorphism."""

    def __init__(self):
        self.buckets = {}
        self.reps = []      # list of (F, L, G, count_labelled)

    def add(self, F, L, n=N):
        G = incidence_graph(F, L, n)
        h = wl_hash(G)
        for idx in self.buckets.get(h, []):
            if isomorphic(self.reps[idx][2], G):
                self.reps[idx][3] += 1
                return idx
        self.reps.append([F, L, G, 1])
        self.buckets.setdefault(h, []).append(len(self.reps) - 1)
        return len(self.reps) - 1

    def classes(self):
        return [(r[0], r[1], r[3]) for r in self.reps]


def automorphisms(F, L=(), n=N):
    """All permutations of P preserving (F, L) (as a list of tuples img[p])."""
    G = incidence_graph(F, L, n)
    nm = isomorphism.categorical_node_match('c', -1)
    gm = isomorphism.GraphMatcher(G, G, node_match=nm)
    out = set()
    for m in gm.isomorphisms_iter():
        out.add(tuple(m[('p', p)][1] for p in range(n)))
    return sorted(out)


# ----------------------------------------------------------------------------------------------
# real linear spaces forced by (F, L)
def derived9(F, p):
    """lines of the derived (n-1)-point set at p (as frozensets of points != p)."""
    return [frozenset(B) - {p} for B in F if p in B]


def derived_Q(F, L, p, inf='oo'):
    """lines (>= 3 points) of the 10-point set Q_p = (P - p) u {oo'}."""
    Lset = set(frozenset(S) for S in L)
    out = []
    for B in F:
        if p in B:
            B = frozenset(B)
            d = B - {p}
            if B in Lset:
                d = d | {inf}
            out.append(d)
    for S in L:
        if len(S) == 3 and p in S:
            out.append((frozenset(S) - {p}) | {inf})
    return out


def line_space(F, L):
    return [frozenset(S) for S in L]


# ----------------------------------------------------------------------------------------------
# filters on a real linear space (points, lines); lines pairwise share <= 1 point
def hereditary_sg_violation(points, lines, min_size=3):
    """Return a subset U (|U| >= min_size) not contained in a single line such that every pair of U
    lies on a line with a third point of U (i.e. U has no ordinary line) -- a Sylvester-Gallai
    violation -- or None.  Complete check over all subsets of size >= min_size."""
    points = sorted(points, key=str)
    lines = [frozenset(l) for l in lines]
    for k in range(min_size, len(points) + 1):
        for U in itertools.combinations(points, k):
            Us = frozenset(U)
            if any(Us <= l for l in lines):
                continue
            ok = True
            for a, b in itertools.combinations(U, 2):
                # pair {a,b} must be on a line with >= 3 points of U to be non-ordinary
                found = False
                for l in lines:
                    if a in l and b in l:
                        if len(l & Us) >= 3:
                            found = True
                        break
                if not found:
                    ok = False
                    break
            if ok:
                return U
    return None


def sg_pair_cap_violation(points, lines):
    """Return True if the lines cover more than C(m,2) - 1 pairs (m = #points) and the points are
    not all on one line."""
    m = len(points)
    lines = [frozenset(l) for l in lines]
    if any(len(l) == m for l in lines):
        return False
    covered = sum(comb(len(l), 2) for l in lines)
    return covered > comb(m, 2) - 1


def mk_violation(points, lines):
    """(8_3) lemma: no 8 distinct real points carry 8 three-point lines.  Return an 8-subset E such
    that >= 8 lines meet E in exactly 3 points, or None."""
    points = sorted(points, key=str)
    lines = [frozenset(l) for l in lines]
    for E in itertools.combinations(points, 8):
        Es = frozenset(E)
        c = sum(1 for l in lines if len(l & Es) == 3)
        if c >= 8:
            return E
    return None


def t3_cap_violation(points, lines, caps):
    """orchard caps: caps = {m: max number of exact-3-point lines on m points}.  Checks every
    subset U of size m in caps (lines restricted to U with exactly 3 points)."""
    points = sorted(points, key=str)
    lines = [frozenset(l) for l in lines]
    for m, cap in caps.items():
        if m > len(points):
            continue
        for U in itertools.combinations(points, m):
            Us = frozenset(U)
            c = sum(1 for l in lines if len(l & Us) == 3)
            if c > cap:
                return (m, U, c)
    return None


PROVED_T3 = {7: 6, 8: 7, 9: 10, 10: 13}   # 7: Fano not real (SG); 8: (8_3); 9: STS(9)-line; 10: parity
CITED_T3 = {7: 6, 8: 7, 9: 10, 10: 12}    # Burr-Gruenbaum-Sloane 1974


def all_spaces(F, L, n=N):
    """yield (name, points, lines) for the 2n+1 real linear spaces forced by (F, L)."""
    for p in range(n):
        pts = [q for q in range(n) if q != p]
        yield ('d9@%d' % p, pts, derived9(F, p))
    for p in range(n):
        pts = [q for q in range(n) if q != p] + ['oo']
        yield ('Q@%d' % p, pts, derived_Q(F, L, p))
    yield ('lines', list(range(n)), line_space(F, L))


def run_filters(F, L, n=N, cited=False):
    """Return list of (space, filter, witness) violations."""
    out = []
    for name, pts, lines in all_spaces(F, L, n):
        if sg_pair_cap_violation(pts, lines):
            out.append((name, 'SGcap', None))
        w = hereditary_sg_violation(pts, lines)
        if w is not None:
            out.append((name, 'HSG', w))
        w = mk_violation(pts, lines)
        if w is not None:
            out.append((name, 'MK', w))
        w = t3_cap_violation(pts, lines, CITED_T3 if cited else PROVED_T3)
        if w is not None:
            out.append((name, 'T3', w))
    w = miquel_violation(F, L, n)
    if w is not None:
        out.append(('mobius', 'MIQ', w))
    return out


# ----------------------------------------------------------------------------------------------
# Miquel closure on P u {oo}
def mobius_blocks(F, L, n=N):
    """blocks of the Möbius structure on P u {oo}: rich blocks of F (as circles) and lines with oo
    appended.  Only blocks of size >= 4 in P u {oo} are returned as 'rich' (F members and lines)."""
    Lset = set(frozenset(S) for S in L)
    out = []
    for B in F:
        B = frozenset(B)
        if B in Lset:
            out.append(B | {'oo'})
        else:
            out.append(B)
    for S in L:
        S = frozenset(S)
        if len(S) == 3:
            out.append(S | {'oo'})
    return out


CUBE_FACES = [((0, 1, 2, 3), (4, 5, 6, 7)), ((0, 1, 5, 4), (3, 2, 6, 7)), ((0, 3, 7, 4), (1, 2, 6, 5))]
# vertex i of the cube: bits (x,y,z) -> 0:(000) 1:(100) 2:(110) 3:(010) 4:(001) 5:(101) 6:(111) 7:(011)


def miquel_violation(F, L, n=N):
    """Miquel: 8 distinct points of the Möbius plane labelled by a cube cannot have exactly 5
    concyclic faces.  Search for a cube (8 points of P u {oo}) with 5 faces contained in blocks and the
    6th face not contained in any block.  Return the offending labelling or None.
    Implementation: choose two disjoint 4-subsets A, B contained in blocks (the 'top' and 'bottom'
    faces) and a bijection between them (4 side faces = {a_i, a_{i+1}, b_{i+1}, b_i}); the cyclic
    orders of A and B must be consistent; enumerate all cyclic orders of A (3 up to reversal) and all
    bijections A -> B (24)."""
    blocks = mobius_blocks(F, L, n)
    fours = set()
    for Bk in blocks:
        for Q in itertools.combinations(sorted(Bk, key=str), 4):
            fours.add(frozenset(Q))
    fl = sorted(fours, key=lambda s: sorted(map(str, s)))
    for i in range(len(fl)):
        A = fl[i]
        for j in range(len(fl)):
            if i == j:
                continue
            B = fl[j]
            if A & B:
                continue
            a = sorted(A, key=str)
            # cyclic orders of A up to rotation/reflection: fix a[0] first, then 3 orders
            for perm in itertools.permutations(a[1:]):
                if str(perm[0]) > str(perm[2]):
                    continue
                cyc = (a[0],) + perm
                for bperm in itertools.permutations(sorted(B, key=str)):
                    faces = [frozenset(cyc), frozenset(bperm)]
                    for k in range(4):
                        faces.append(frozenset((cyc[k], cyc[(k + 1) % 4], bperm[(k + 1) % 4], bperm[k])))
                    inb = [f in fours for f in faces]
                    if sum(inb) == 5:
                        return (cyc, bperm, [sorted(f, key=str) for f in faces], inb)
    return None
