"""Erdős #506 — independent CP-SAT model of the Möbius combinatorial relaxation (n10-cpsat/v2 agent, 2026-09-05).
Written from scratch; nothing imported from the other experiment directories.

Setting.  P = {0..n-1} is a finite set in the real Möbius plane (plane + infinity), not all on one circle/line.
A *rich block* is a circle-or-line through >= 4 points of P (as a subset of P); a *line* is a block through the
point at infinity, i.e. a collinear subset of size >= 3.  circles(P) = C(n,3) - D - l with
D = sum_{rich B} (C(|B|,3) - 1) and l = number of lines.  We look for (F, L) with D + l >= threshold.

Necessary conditions encoded (all are theorems unless flagged):
 (T)   every triple lies in <= 1 rich block; a 3-line is a triple in no rich block.
 (L)   a line of size >= 4 is a rich block; two lines share <= 1 point.
 (PD)  the blocks through a pair {p,q} partition the other n-2 points: sum_{B ⊇ {p,q}} (|B|-2) x_B <= n-2;
       the lines through p partition the other points: sum_{S ∋ p} (|S|-1) y_S <= n-1.
 (SG)  derived Sylvester–Gallai: inverting about p ∈ P turns the rich blocks through p into the >= 3-point lines
       of a non-collinear real (n-1)-point set, which has >= o(n-1) ordinary lines; hence the rich blocks through
       p cover <= C(n-1,2) - o(n-1) pairs.  Lines: <= C(n,2) - o(n) pairs.  o(m) >= 1 always [Sylvester–Gallai];
       o(m) >= 3m/7 [Kelly–Moser 1958]; exact table o(9)=6, o(10)=5 [Crowe–McKee 1968] — mode 'table' is CITED.
 (HSG) hereditary SG on r-subsets S of a derived point set (r = 7, 8) and of P (r = 7, 8, 9): the >= 3-point lines
       of S (restrictions of derived lines meeting S in >= 3 points) cover <= C(r,2) - 1 pairs.  (r = 7 excludes
       a derived Fano plane.)
 (MK)  no 8 real points carry 8 three-point lines (the (8_3) Möbius–Kantor configuration is not realisable over R;
       re-verified in facts.py): for every 8-subset S of a derived point set, at most 7 derived lines meet S in
       exactly 3 points; same for the lines of P.  Consequence (facts.py): 9 real points carry <= 10 three-point
       lines, so <= 10 four-blocks through any point (O9), and <= 10 three-lines inside any 9-subset of P.
 (E11) the same at p for the 11-point Möbius set P ∪ {∞}: after inversion about p the lines through p become the
       3-point lines (S - p) ∪ {∞'} of a real n-point set; SG, HSG, MK apply to it (subsets containing ∞').
 (O10) [OPTIONAL, cited: orchard number t3(10) = 12, Burr–Grünbaum–Sloane 1974] <= 12 three-lines in total and, at
       every p, (#four-blocks through p) + (#three-lines through p) <= 12.
Symmetry breaking (valid for enumeration up to isomorphism): a chosen family of big blocks (the 'skeleton') is
fixed; within every cell of points having the same membership pattern in the skeleton the degree key is
non-increasing.
"""
import itertools, math
from ortools.sat.python import cp_model

O_MODES = {'sg': lambda m: 1, 'km': lambda m: -(-3 * m // 7), 'table': lambda m: {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5, 11: 6, 12: 6, 13: 6, 14: 7}[m]}
T3_PROVED = {3: 1, 4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10}   # from packing + SG(Fano) + MK (facts.py)
T3_CITED = {10: 12, 11: 16, 12: 19}                       # Burr–Grünbaum–Sloane 1974


def deficit(k):
    return math.comb(k, 3) - 1


class Relaxation:
    def __init__(self, n=10, sizes=(4, 5, 6), o='sg', threshold=None, skeleton=None, forbid_other_big=True, forbid_sizes=None,
                 hsg=True, mk=True, orchard9=True, e11=True, orchard10=False, pairdeg=True, degree_order=True):
        self.n = n
        PTS = self.PTS = tuple(range(n))
        self.sizes = sizes
        self.o = O_MODES[o]
        self.omode = o
        m = self.m = cp_model.CpModel()
        self.blocks = [frozenset(c) for k in sizes for c in itertools.combinations(PTS, k)]
        self.lines = [frozenset(c) for k in (3,) + tuple(s for s in sizes) for c in itertools.combinations(PTS, k)]
        self.x = {B: m.NewBoolVar('x' + ''.join(map(str, sorted(B)))) for B in self.blocks}
        self.y = {S: m.NewBoolVar('y' + ''.join(map(str, sorted(S)))) for S in self.lines}
        x, y = self.x, self.y
        self.bl_by_point = {p: [B for B in self.blocks if p in B] for p in PTS}
        self.ln_by_point = {p: [S for S in self.lines if p in S] for p in PTS}
        self.ncons = {}
        def cnt(tag):
            self.ncons[tag] = self.ncons.get(tag, 0) + 1
        # (T)
        for T in itertools.combinations(PTS, 3):
            T = frozenset(T)
            m.AddAtMostOne([y[T]] + [x[B] for B in self.blocks if T <= B]); cnt('T')
        # (L)
        for S in self.lines:
            if len(S) >= 4:
                m.AddImplication(y[S], x[S]); cnt('L1')
        for Q in itertools.combinations(PTS, 2):
            Q = frozenset(Q)
            m.AddAtMostOne([y[S] for S in self.lines if Q <= S]); cnt('L2')
            if pairdeg:
                m.Add(sum((len(B) - 2) * x[B] for B in self.blocks if Q <= B) <= n - 2); cnt('PD')
        if pairdeg:
            for p in PTS:
                m.Add(sum((len(S) - 1) * y[S] for S in self.ln_by_point[p]) <= n - 1); cnt('PDl')
        # (LPK) local packing at a point (facts.py F1 and the extra table: with (d6, d5) derived 5-/4-lines on the
        # 9 other points the number of derived 3-lines is <= 10 (MK), 8, 7, 4 for d5 = 0..3 (d6 = 0), <= 6, 3 for
        # (d6, d5) = (1, 0), (1, 1), and 0 for (2, 0); (1, 2), (2, 1) impossible).  Valid linear cuts (n = 10 only):
        if n == 10:
            for p in PTS:
                d4 = sum(x[B] for B in self.bl_by_point[p] if len(B) == 4)
                d5 = sum(x[B] for B in self.bl_by_point[p] if len(B) == 5)
                d6 = sum(x[B] for B in self.bl_by_point[p] if len(B) == 6)
                m.Add(2 * d4 + 3 * d5 + 8 * d6 <= 20); cnt('LPK')
                m.Add(d5 + 2 * d6 <= 3); cnt('LPK')
        # (SG)
        for p in PTS:
            m.Add(sum(math.comb(len(B) - 1, 2) * x[B] for B in self.bl_by_point[p]) <= math.comb(n - 1, 2) - self.o(n - 1)); cnt('SG')
        m.Add(sum(math.comb(len(S), 2) * y[S] for S in self.lines) <= math.comb(n, 2) - self.o(n)); cnt('SGl')
        # (HSG) + (MK) for derived structures and for the lines
        def add_subset_caps(point_set, linevars, tag):
            """linevars: list of (set_of_points, var); hereditary SG (r=7,8) and MK (r=8), t3(9)<=10 (r=9)."""
            for r in (7, 8, 9):
                if r > len(point_set):
                    continue
                for S in itertools.combinations(sorted(point_set), r):
                    S = frozenset(S)
                    terms3 = [v for (Ls, v) in linevars if len(Ls & S) == 3]
                    if hsg and r in (7, 8):
                        terms = [(math.comb(len(Ls & S), 2), v) for (Ls, v) in linevars if len(Ls & S) >= 3]
                        if sum(c for c, _ in terms) > math.comb(r, 2) - 1:
                            m.Add(sum(c * v for c, v in terms) <= math.comb(r, 2) - 1); cnt(tag + f'HSG{r}')
                    if mk and r == 8 and len(terms3) > 7:
                        m.Add(sum(terms3) <= 7); cnt(tag + 'MK8')
                    if mk and orchard9 and r == 9 and len(terms3) > 10:
                        m.Add(sum(terms3) <= 10); cnt(tag + 'O9')
        for p in PTS:
            others = frozenset(q for q in PTS if q != p)
            dl = [(B - {p}, x[B]) for B in self.bl_by_point[p]]
            add_subset_caps(others, dl, 'd')
            if e11:
                # derived structure at p of the 11-point Möbius set P ∪ {∞}: a real n-point set (P - p) ∪ {∞'}.
                # Its >= 3-point lines: for a rich block B ∋ p the line B - p, extended by ∞' iff B is also a line
                # (y_B = 1, which implies x_B = 1); for a 3-line S ∋ p the 3-point line (S - p) ∪ {∞'}.
                # (A rich line through p is ONE line of the derived set — no double counting.)
                INF = n
                rich = self.bl_by_point[p]
                three = [S for S in self.ln_by_point[p] if len(S) == 3]
                m.Add(sum(math.comb(len(B) - 1, 2) * x[B] + (len(B) - 1) * y[B] for B in rich)
                      + sum(3 * y[S] for S in three) <= math.comb(n, 2) - self.o(n)); cnt('E11SG')
                for r in (7, 8, 9):
                    for S0 in itertools.combinations(sorted(others), r - 1):
                        S0 = frozenset(S0)
                        pair_terms = []; three_terms = []
                        for B in rich:
                            k = len((B - {p}) & S0)
                            if k >= 3:
                                pair_terms += [(math.comb(k, 2), x[B]), (k, y[B])]
                                if k == 3:
                                    three_terms += [(1, x[B]), (-1, y[B])]
                            elif k == 2:
                                pair_terms.append((3, y[B])); three_terms.append((1, y[B]))
                        for S in three:
                            if len((S - {p}) & S0) == 2:
                                pair_terms.append((3, y[S])); three_terms.append((1, y[S]))
                        if hsg and r in (7, 8) and sum(c for c, _ in pair_terms if c > 0) > math.comb(r, 2) - 1:
                            m.Add(sum(c * v for c, v in pair_terms) <= math.comb(r, 2) - 1); cnt(f'E11HSG{r}')
                        if mk and r == 8 and sum(c for c, _ in three_terms if c > 0) > 7:
                            m.Add(sum(c * v for c, v in three_terms) <= 7); cnt('E11MK8')
                        if mk and orchard9 and r == 9 and sum(c for c, _ in three_terms if c > 0) > 10:
                            m.Add(sum(c * v for c, v in three_terms) <= 10); cnt('E11O9')
                if orchard10:
                    m.Add(sum(x[B] - y[B] for B in rich if len(B) == 4) + sum(y[S] for S in three) <= T3_CITED[10]); cnt('E11O10')
        add_subset_caps(frozenset(PTS), [(S, y[S]) for S in self.lines], 'l')
        if orchard10:
            m.Add(sum(y[S] for S in self.lines if len(S) == 3) <= T3_CITED[10]); cnt('O10')
        # skeleton: fixed big blocks
        self.skeleton = None
        if skeleton is not None:
            self.skeleton = [frozenset(B) for B in skeleton]
            for B in self.skeleton:
                m.Add(x[B] == 1)
            if forbid_other_big:
                fs = set(forbid_sizes) if forbid_sizes is not None else {s for s in sizes if s >= 5}
                for B in self.blocks:
                    if len(B) in fs and B not in self.skeleton:
                        m.Add(x[B] == 0)
            # cells and degree ordering
            if degree_order:
                cells = {}
                for p in PTS:
                    key = tuple(p in B for B in self.skeleton)
                    cells.setdefault(key, []).append(p)
                self.cells = list(cells.values())
                for cell in self.cells:
                    for a, b in zip(cell, cell[1:]):
                        m.Add(self.degkey(a) >= self.degkey(b)); cnt('SYM')
        self.D = sum(deficit(len(B)) * x[B] for B in self.blocks)
        self.l = sum(y[S] for S in self.lines)
        self.obj = self.D + self.l
        if threshold is not None:
            m.Add(self.obj >= threshold)

    def degkey(self, p):
        # blocks of size s through p weighted by 64^(s-4) (larger blocks dominate), plus lines through p in the lowest digit
        return sum(64 ** (len(B) - 4 + 1) * self.x[B] for B in self.bl_by_point[p]) + sum(self.y[S] for S in self.ln_by_point[p])

    def extract(self, val):
        F = [B for B in self.blocks if val(self.x[B])]
        L = [S for S in self.lines if val(self.y[S])]
        return F, L


def count(n, F, L):
    return math.comb(n, 3) - sum(deficit(len(B)) for B in F) - len(L)


# ------------------------------------------------------------------ isomorphism bookkeeping (own code, networkx used only for VF2)
import networkx as nx
from collections import Counter

def incidence_graph(n, F, L):
    G = nx.Graph()
    for p in range(n):
        G.add_node(('p', p), c=0)
    for i, B in enumerate(F):
        G.add_node(('b', i), c=10 + len(B))
        for p in B:
            G.add_edge(('b', i), ('p', p))
    for i, S in enumerate(L):
        G.add_node(('l', i), c=30 + len(S))
        for p in S:
            G.add_edge(('l', i), ('p', p))
    return G


def invariant(n, F, L):
    prof = []
    for p in range(n):
        c = Counter()
        for B in F:
            if p in B: c[('b', len(B))] += 1
        for S in L:
            if p in S: c[('l', len(S))] += 1
        prof.append(tuple(sorted(c.items())))
    # pair profile: multiset over pairs of (size of common block, in a common line?)
    pairs = Counter()
    for a, b in itertools.combinations(range(n), 2):
        bs = tuple(sorted(len(B) for B in F if a in B and b in B))
        ls = tuple(sorted(len(S) for S in L if a in S and b in S))
        pairs[(bs, ls)] += 1
    return (tuple(sorted(len(B) for B in F)), tuple(sorted(len(S) for S in L)), tuple(sorted(prof)), tuple(sorted(pairs.items())),
            nx.weisfeiler_lehman_graph_hash(incidence_graph(n, F, L), node_attr='c', iterations=4))


NM = nx.algorithms.isomorphism.categorical_node_match('c', 0)

class ClassStore:
    def __init__(self, n=10):
        self.n = n; self.buckets = {}; self.classes = []; self.total = 0
    def add(self, F, L):
        F = [frozenset(B) for B in F]; L = [frozenset(S) for S in L]
        self.total += 1
        inv = invariant(self.n, F, L)
        G = incidence_graph(self.n, F, L)
        for idx in self.buckets.get(inv, []):
            if nx.is_isomorphic(self.classes[idx]['G'], G, node_match=NM):
                self.classes[idx]['copies'] += 1
                return idx, False
        self.buckets.setdefault(inv, []).append(len(self.classes))
        self.classes.append({'F': F, 'L': L, 'inv': inv, 'copies': 1, 'G': G})
        return len(self.classes) - 1, True
    def records(self):
        return [{'blocks': sorted(sorted(B) for B in c['F']), 'lines': sorted(sorted(S) for S in c['L']),
                 'copies': c['copies'], 'count': count(self.n, c['F'], c['L'])} for c in self.classes]


def automorphisms(n, F, L=()):
    """all permutations of the points preserving (F, L), as tuples."""
    G = incidence_graph(n, F, list(L))
    GM = nx.algorithms.isomorphism.GraphMatcher(G, G, node_match=NM)
    auts = []
    for iso in GM.isomorphisms_iter():
        auts.append(tuple(iso[('p', p)][1] for p in range(n)))
    return auts
