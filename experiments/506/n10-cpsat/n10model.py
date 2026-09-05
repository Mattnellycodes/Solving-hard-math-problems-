"""Erdős #506, n = 10: CP-SAT model of the Möbius combinatorial relaxation, restricted to
largest block size 5 or 6.  Written from scratch (n10-cpsat agent, 2026-09-05).

Setting.  P = {0..9}.  A *rich block* is a subset of size >= 4 (circle or line through >= 4 points);
a *line* is a subset of size >= 3 that is collinear (block through infinity).  Facts used:
  (T)  every triple lies in <= 1 rich block (two rich blocks share <= 2 points);
  (L)  a line of size >= 4 is a rich block; a 3-line is a triple in no rich block; two lines share
       <= 1 point;
  (SG) derived Sylvester–Gallai at p: sum_{B ∋ p} C(|B|-1,2) x_B <= C(9,2) - o(9)   (o = min number of
       ordinary lines of a non-collinear m-point set; o(m) >= 1 always [Sylvester–Gallai];
       o(9) >= 5, o(10) >= 5 [Csima–Sawyer 6m/13, m != 7]; exact o(9) = 6, o(10) = 5 [Crowe–McKee table]);
  (SGL) line cap: sum_S C(|S|,2) y_S <= C(10,2) - o(10);
  (PD) pair degree: blocks through {p,q} partition the other 8 points: sum_{B ⊇ {p,q}} (|B|-2) x_B <= 8,
       and lines through p: sum_{S ∋ p} (|S|-1) y_S <= 9   (pure packing, valid);
  (F7) no derived Fano plane: for every p and 7-subset S of P\{p}, at most 6 blocks through p meet S
       in exactly 3 points (7 such would be 7 lines on 7 points pairwise sharing <= 1 point, i.e. a Fano
       plane covering all 21 pairs, contradicting SG on S); same for the lines and 7-subsets of P;
  (M8) no derived (8_3): for every p and 8-subset S of P\{p}, at most 7 blocks through p meet S in
       exactly 3 points (8 such are 8 triples on 8 points pairwise sharing <= 1 point = Möbius–Kantor,
       not realisable over R: verified in ../verify_independent/mk_unrealisable2.py and
       ../theory/eight_three_light.py); same for the lines and 8-subsets of P;
  (O9) <= 10 four-blocks through any point.  This is the orchard value t3(9) = 10 of Burr–Grünbaum–
       Sloane (1974); it is ALSO implied by (M8): 11 pairwise <=1-intersecting triples on 9 points
       form AG(2,3) minus a line, which contains an (8_3) on the 8 points off a removed-line point.
       (Kept as an explicit, separately switchable constraint.)
  (O10) <= 12 three-point lines in total: orchard value t3(10) = 12 (Burr–Grünbaum–Sloane 1974),
       CITED, not proved here; switchable, off by default.
  (D5) each point lies on <= 4 five-blocks and <= 2 six-blocks (packing: k blocks of size s through p
       pairwise sharing <= 1 further point need 1 + k(s-1) - ... points; direct check), redundant but
       helps the LP.
Objective/threshold:  D + l >= 88  where D = sum_B (C(|B|,3)-1) x_B, l = number of lines;
circles = 120 - D - l, so D + l >= 88  <=>  circles <= 32.
Cases:  'A': a 6-block exists, fixed to {0..5}; block sizes 4,5,6.
        'B': largest block 5, fixed to {0..4}; block sizes 4,5.
        'C': largest block 4, fixed to {0..3}; block size 4 only (outside the assigned task; extra check).
Symmetry breaking (valid for enumeration up to isomorphism): the stabiliser of the fixed block is
S_k x S_{10-k}; we require a point key (d6, d5, d4) to be non-increasing inside the fixed block and
inside its complement.  Every isomorphism class therefore still has >= 1 labelled representative.
"""
import itertools, math
from collections import Counter
from ortools.sat.python import cp_model

N = 10
PTS = tuple(range(N))
TRIPLES = [frozenset(t) for t in itertools.combinations(PTS, 3)]
PAIRS = [frozenset(t) for t in itertools.combinations(PTS, 2)]

O_EXACT = {9: 6, 10: 5}       # Crowe–McKee table (o(3..14) = 3,3,4,3,3,4,6,5,6,6,6,7)
O_CS = {9: 5, 10: 5}          # Csima–Sawyer 6m/13 rounded up (m != 7)
O_SG = {9: 1, 10: 1}          # Sylvester–Gallai only
O_MODES = {'table': O_EXACT, 'cs': O_CS, 'sg': O_SG}


def deficit(k):
    return math.comb(k, 3) - 1


class Model:
    def __init__(self, case, omode='table', fano=True, mk=True, orchard9=True, orchard10=False,
                 pairdeg=True, degree_order=True, threshold=88, extra_fixed=None, bundle=False, eleven=False):
        assert case in ('A', 'B', 'C')
        self.case = case
        self.o = O_MODES[omode]
        self.omode = omode
        sizes = {'A': (4, 5, 6), 'B': (4, 5), 'C': (4,)}[case]
        lsizes = (3,) + sizes
        self.sizes, self.lsizes = sizes, lsizes
        m = cp_model.CpModel()
        self.m = m
        self.blocks = [frozenset(c) for k in sizes for c in itertools.combinations(PTS, k)]
        self.lines = [frozenset(c) for k in lsizes for c in itertools.combinations(PTS, k)]
        self.x = {B: m.NewBoolVar('x' + ''.join(map(str, sorted(B)))) for B in self.blocks}
        self.y = {S: m.NewBoolVar('y' + ''.join(map(str, sorted(S)))) for S in self.lines}
        x, y = self.x, self.y
        # index: blocks containing a given triple / pair / point
        self.bl_by_triple = {T: [B for B in self.blocks if T <= B] for T in TRIPLES}
        self.bl_by_point = {p: [B for B in self.blocks if p in B] for p in PTS}
        self.ln_by_pair = {Q: [S for S in self.lines if Q <= S] for Q in PAIRS}
        self.ln_by_point = {p: [S for S in self.lines if p in S] for p in PTS}
        self.bl_by_pair = {Q: [B for B in self.blocks if Q <= B] for Q in PAIRS}
        # (T) + 3-lines are uncovered triples
        for T in TRIPLES:
            m.AddAtMostOne([y[T]] + [x[B] for B in self.bl_by_triple[T]])
        # (L) lines of size >= 4 are rich blocks; lines pairwise share <= 1 point
        for S in self.lines:
            if len(S) >= 4:
                m.AddImplication(y[S], x[S])
        for Q in PAIRS:
            m.AddAtMostOne([y[S] for S in self.ln_by_pair[Q]])
        # (SG) derived caps
        for p in PTS:
            m.Add(sum(math.comb(len(B) - 1, 2) * x[B] for B in self.bl_by_point[p]) <= math.comb(9, 2) - self.o[9])
        # (SGL)
        m.Add(sum(math.comb(len(S), 2) * y[S] for S in self.lines) <= math.comb(10, 2) - self.o[10])
        # (PD)
        if pairdeg:
            for Q in PAIRS:
                m.Add(sum((len(B) - 2) * x[B] for B in self.bl_by_pair[Q]) <= 8)
            for p in PTS:
                m.Add(sum((len(S) - 1) * y[S] for S in self.ln_by_point[p]) <= 9)
        # (D5)/(DD): local packing bounds at a point, verified by brute force in check_local_bounds.py:
        # with k derived 4-lines on the 9 other points, the number of derived 3-lines is <= 12, 8, 7, 4
        # for k = 0, 1, 2, 3 and k = 4 is impossible.  (k = 0: 10 by (M8)/(O9).)
        for p in PTS:
            d5 = sum(x[B] for B in self.bl_by_point[p] if len(B) == 5)
            d4 = sum(x[B] for B in self.bl_by_point[p] if len(B) == 4)
            if 5 in sizes:
                m.Add(d5 <= 3)
                m.Add(d4 + 3 * d5 <= 13)            # (2,7), (3,4)
                for B in self.bl_by_point[p]:
                    if len(B) == 5:
                        m.Add(d4 + 2 * x[B] <= 10)  # (1,8)
            if 6 in sizes:
                m.Add(sum(x[B] for B in self.bl_by_point[p] if len(B) == 6) <= 2)
        # (F7), (M8) for derived structures and for the lines
        self.n_sub = 0
        for r, cap, flag in ((7, 6, fano), (8, 7, mk)):
            if not flag:
                continue
            for p in PTS:
                others = [q for q in PTS if q != p]
                for S in itertools.combinations(others, r):
                    S = frozenset(S)
                    terms = [x[B] for B in self.bl_by_point[p] if len(B & S) == 3]
                    if len(terms) > cap:
                        m.Add(sum(terms) <= cap); self.n_sub += 1
            for S in itertools.combinations(PTS, r):
                S = frozenset(S)
                terms = [y[L] for L in self.lines if len(L & S) == 3]
                if len(terms) > cap:
                    m.Add(sum(terms) <= cap); self.n_sub += 1
        # (O9), (O10)
        if orchard9:
            for p in PTS:
                m.Add(sum(x[B] for B in self.bl_by_point[p] if len(B) == 4) <= 10)
        if orchard10:
            m.Add(sum(y[S] for S in self.lines if len(S) == 3) <= 12)
        # (BT) bundle theorem of the real Möbius plane (sphere model): for four pairwise disjoint pairs
        # P1..P4 of points of P ∪ {∞}, if five of the six quadruples Pi ∪ Pj lie in blocks, so does the
        # sixth.  Proof for the pairs {a,b},{c,d} ⊂ Π-block and {r,r'},{s,s'}: X = ab ∩ cd lies in the
        # planes (abrr'), (cdrr'), hence on line rr'; likewise on ss'; so rr', ss' meet and r,r',s,s' are
        # coplanar, i.e. concyclic.  By symmetry any 5 imply the 6th.  Quadruples through ∞ are lines.
        self.n_bundle = 0
        if bundle:
            INF = N
            allpts = list(PTS) + [INF]
            quad_vars = {}
            def qv(Q):
                Q = frozenset(Q)
                if Q not in quad_vars:
                    if INF in Q:
                        T = Q - {INF}
                        quad_vars[Q] = [y[S] for S in self.lines if T <= S]
                    else:
                        quad_vars[Q] = [x[B] for B in self.blocks if Q <= B]
                return quad_vars[Q]
            pairs11 = [frozenset(t) for t in itertools.combinations(allpts, 2)]
            def rec(chosen, start, used):
                if len(chosen) == 4:
                    terms = []
                    for Pi, Pj in itertools.combinations(chosen, 2):
                        terms += qv(Pi | Pj)
                    z = m.NewBoolVar('')
                    s = sum(terms)
                    m.Add(s <= 4 + 2 * z)
                    m.Add(s >= 6 * z)
                    self.n_bundle += 1
                    return
                for i in range(start, len(pairs11)):
                    Pi = pairs11[i]
                    if Pi & used:
                        continue
                    rec(chosen + [Pi], i + 1, used | Pi)
            rec([], 0, frozenset())
        # (E11) derived caps in the 11-point Möbius set P ∪ {∞}: the derived structure at p is a real 10-point
        # set whose lines are B - p (rich B ∋ p) and (S - p) ∪ {∞'} (lines S ∋ p): pairs <= 45 - o(10),
        # <= 12 three-point lines (t3(10) <= 12, proved here), <= 6 / <= 7 three-lines meeting a 7-/8-subset
        # containing ∞' in exactly 3 points (Fano / (8_3)).
        self.n_e11 = 0
        if eleven:
            for p in PTS:
                m.Add(sum(math.comb(len(B) - 1, 2) * x[B] for B in self.bl_by_point[p])
                      + sum(math.comb(len(S), 2) * y[S] for S in self.ln_by_point[p]) <= math.comb(10, 2) - self.o[10])
                m.Add(sum(x[B] for B in self.bl_by_point[p] if len(B) == 4)
                      + sum(y[S] for S in self.ln_by_point[p] if len(S) == 3) <= 12)
                self.n_e11 += 2
                others = [q for q in PTS if q != p]
                for r, cap in ((7, 6), (8, 7)):
                    for S in itertools.combinations(others, r - 1):      # subsets containing ∞'
                        S = frozenset(S)
                        terms = [x[B] for B in self.bl_by_point[p] if len(B & S) == 3] + \
                                [y[T] for T in self.ln_by_point[p] if len(T & S) == 2]
                        if len(terms) > cap:
                            m.Add(sum(terms) <= cap); self.n_e11 += 1
        # fixed largest block
        k = {'A': 6, 'B': 5, 'C': 4}[case]
        self.fixed = frozenset(range(k))
        m.Add(x[self.fixed] == 1)
        if extra_fixed:
            for B, v in extra_fixed.items():
                m.Add(x[frozenset(B)] == v)
        # degree keys and symmetry breaking
        self.d = {}
        for p in PTS:
            for s in sizes:
                self.d[p, s] = sum(x[B] for B in self.bl_by_point[p] if len(B) == s)
        def key(p):
            return sum(self.d[p, s] * (100 ** (s - 4)) for s in sizes)
        if degree_order:
            for grp in (list(range(k)), list(range(k, N))):
                for a, b in zip(grp, grp[1:]):
                    m.Add(key(a) >= key(b))
        # objective
        self.D = sum(deficit(len(B)) * x[B] for B in self.blocks)
        self.l = sum(y[S] for S in self.lines)
        self.obj = self.D + self.l
        if threshold is not None:
            m.Add(self.obj >= threshold)

    def extract(self, val):
        F = [B for B in self.blocks if val(self.x[B])]
        L = [S for S in self.lines if val(self.y[S])]
        return F, L


def structure_count(F, L):
    return math.comb(N, 3) - sum(deficit(len(B)) for B in F) - len(L)


# ---------------------------------------------------------------------------------------------
# isomorphism-class bookkeeping (networkx); the invariant is only a bucket key, the isomorphism
# test is exact.
import networkx as nx

def incidence_graph(F, L):
    G = nx.Graph()
    for p in PTS:
        G.add_node(('p', p), c=0)
    for i, B in enumerate(F):
        G.add_node(('b', i), c=10 + len(B))
        for p in B:
            G.add_edge(('b', i), ('p', p))
    for i, S in enumerate(L):
        G.add_node(('l', i), c=20 + len(S))
        for p in S:
            G.add_edge(('l', i), ('p', p))
    return G


def invariant(F, L):
    deg = {p: Counter() for p in PTS}
    for B in F:
        for p in B:
            deg[p][('b', len(B))] += 1
    for S in L:
        for p in S:
            deg[p][('l', len(S))] += 1
    prof = tuple(sorted(tuple(sorted(deg[p].items())) for p in PTS))
    return (tuple(sorted(len(B) for B in F)), tuple(sorted(len(S) for S in L)), prof,
            nx.weisfeiler_lehman_graph_hash(incidence_graph(F, L), node_attr='c', iterations=4))


class ClassStore:
    """Collects labelled structures and keeps one representative per isomorphism class."""
    def __init__(self):
        self.buckets = {}
        self.classes = []      # list of (F, L, inv, count_of_labelled_copies)
        self.total = 0

    def add(self, F, L):
        self.total += 1
        inv = invariant(F, L)
        G = incidence_graph(F, L)
        nm = nx.algorithms.isomorphism.categorical_node_match('c', 0)
        for idx in self.buckets.get(inv, []):
            if nx.is_isomorphic(self.classes[idx][4], G, node_match=nm):
                self.classes[idx][3] += 1
                return False
        self.buckets.setdefault(inv, []).append(len(self.classes))
        self.classes.append([F, L, inv, 1, G])
        return True
