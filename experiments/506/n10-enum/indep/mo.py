"""mo.py -- independent two-phase orderly enumeration of abstract Moebius block structures
(Erdos #506, lower bound for c(10)).  Written from scratch (n10-enum/indep agent, 2026-09-05);
shares no code with ../orderly2.py, ../../theory/mobius_enum*.py or the CP-SAT models.

Objects.  Points 0..n-1 (bitmasks).  A *rich block* is a subset of size 4..n-1; a family F of rich
blocks must be pairwise <=2-intersecting.  D(F) = sum (C(|B|,3)-1).  *Lines* are blocks of F or
uncovered triples, pairwise <=1-intersecting.  circles = C(n,3) - D(F) - |L|.

Necessary conditions used (per mode):
  'none'  : intersection conditions only (validation runs).
  'sg'    : derived Sylvester-Gallai with o(m) >= 1 only:  for every point p,
            sum_{B in F, p in B} C(|B|-1,2) <= C(n-1,2) - 1   (the derived rich lines at p cover
            pairwise disjoint pair sets and miss at least one ordinary pair);
            lines: sum C(|L|,2) <= C(n,2) - 1.
  'table' : the same with the known values o(m) (min number of ordinary lines, cited table) and the
            orchard numbers t3(m) (Burr-Gruenbaum-Sloane 1974, cited): at most t3(n-1) four-blocks
            through a point, at most t3(n) three-point lines.

Two-phase orderly generation.  Phase 1: one representative per isomorphism class of families of
*big* blocks (size >= 5) (level-by-level extension + explicit isomorphism rejection).  Phase 2: for a
fixed representative R, families R + F4 (F4 = set of 4-blocks) are isomorphic iff F4, F4' are in the
same orbit of Aut(R); we generate exactly the F4 that are lexicographically minimal (as sorted index
lists) in their Aut(R)-orbit, adding 4-blocks in increasing index order.  Soundness (Read-Faradzhev):
if F4 is orbit-minimal then F4 minus its largest element is orbit-minimal (see NOTES in the report);
hence DFS over children with larger indices + minimality test is exhaustive and duplicate-free.
"""
import itertools, math, sys, time, json
import numpy as np

popcount = int.bit_count
INF = 10 ** 9

# minimum number of ordinary lines o(m) (Kelly-Moser / Csima-Sawyer / known small cases; cited) and
# orchard numbers t3(m) = max number of 3-point lines (Burr-Gruenbaum-Sloane 1974; cited)
O_MIN = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5, 11: 6, 12: 6, 13: 6, 14: 7}
T3 = {3: 1, 4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10, 10: 12, 11: 16, 12: 19}


def masks_of_size(n, k):
    return sorted(sum(1 << i for i in c) for c in itertools.combinations(range(n), k))


def bits(m):
    out = []
    while m:
        low = m & -m
        out.append(low.bit_length() - 1)
        m ^= low
    return out


def mask(pts):
    return sum(1 << p for p in pts)


def apply_perm(sigma, m):
    r = 0
    while m:
        low = m & -m
        r |= 1 << sigma[low.bit_length() - 1]
        m ^= low
    return r


def deficit(F):
    return sum(math.comb(popcount(b), 3) - 1 for b in F)


class Caps:
    def __init__(self, n, mode):
        self.n, self.mode = n, mode
        if mode == 'none':
            self.capD = self.capL = self.max4 = self.t3n = self.ellmax = INF
            return
        if mode == 'sg':
            o = lambda m: 1
        elif mode == 'table':
            o = lambda m: O_MIN[m]
        else:
            raise ValueError(mode)
        self.capD = math.comb(n - 1, 2) - o(n - 1)      # pair budget of the derived structure at a point
        self.capL = math.comb(n, 2) - o(n)              # pair budget of the lines
        self.max4 = self.capD // 3                      # 4-blocks through a point
        self.t3n = INF                                  # 3-point lines at infinity
        if mode == 'table':
            self.max4 = min(self.max4, T3[n - 1])
            self.t3n = T3[n]
        k3 = min(self.t3n, self.capL // 3)
        self.ellmax = k3 + (self.capL - 3 * k3) // 6    # global upper bound for the number of lines

    def __repr__(self):
        return (f"Caps(n={self.n}, mode={self.mode}, capD={self.capD}, max4={self.max4}, "
                f"capL={self.capL}, t3n={self.t3n}, ellmax={self.ellmax})")


def point_usage(n, F):
    return [sum(math.comb(popcount(b) - 1, 2) for b in F if b >> p & 1) for p in range(n)]


def usage_ok(n, F, caps):
    return max(point_usage(n, F)) <= caps.capD if F else True


# ----------------------------------------------------------------------------------------------
# isomorphisms / automorphisms of block families (backtracking on point images)
# ----------------------------------------------------------------------------------------------
def _profiles(n, F):
    prof = [tuple(sorted(popcount(b) for b in F if b >> p & 1)) for p in range(n)]
    pair = [[None] * n for _ in range(n)]
    for p in range(n):
        for q in range(n):
            if p != q:
                pair[p][q] = tuple(sorted(popcount(b) for b in F if (b >> p & 1) and (b >> q & 1)))
    return prof, pair


def isomorphisms(n, F, G):
    """Generate every permutation sigma (tuple, sigma[p] = image of p) with sigma(F) == G."""
    F = sorted(set(F)); G = sorted(set(G)); Gset = set(G)
    if len(F) != len(G) or sorted(map(popcount, F)) != sorted(map(popcount, G)):
        return
    profF, pairF = _profiles(n, F); profG, pairG = _profiles(n, G)
    if sorted(profF) != sorted(profG):
        return
    order = sorted(range(n), key=lambda p: (-len(profF[p]), p))
    sigma = [-1] * n; used = [False] * n
    Gsizes = {}
    for g in G:
        Gsizes.setdefault(popcount(g), []).append(g)

    def partial_ok(i):
        assigned = 0
        for j in range(i + 1):
            assigned |= 1 << order[j]
        for b in F:
            part = b & assigned
            if not part:
                continue
            img = apply_perm(sigma, part)
            if not any(img & g == img for g in Gsizes[popcount(b)]):
                return False
        return True

    def rec(i):
        if i == n:
            if {apply_perm(sigma, b) for b in F} == Gset:
                yield tuple(sigma)
            return
        p = order[i]
        for q in range(n):
            if used[q] or profG[q] != profF[p]:
                continue
            ok = True
            for j in range(i):
                p2 = order[j]
                if pairF[p][p2] != pairG[q][sigma[p2]]:
                    ok = False; break
            if not ok:
                continue
            sigma[p] = q; used[q] = True
            if partial_ok(i):
                yield from rec(i + 1)
            sigma[p] = -1; used[q] = False
    yield from rec(0)


def is_isomorphic(n, F, G):
    return next(isomorphisms(n, F, G), None) is not None


def automorphisms(n, F):
    return list(isomorphisms(n, F, F))


def invariant(n, F):
    prof = sorted(tuple(sorted(popcount(b) for b in F if b >> p & 1)) for p in range(n))
    inter = sorted(popcount(a & b) for a, b in itertools.combinations(F, 2))
    return (tuple(prof), tuple(inter))


# ----------------------------------------------------------------------------------------------
# phase 1: representatives of families of big blocks
# ----------------------------------------------------------------------------------------------
def phase1(n, caps, sizes, verbose=False):
    big = [m for s in sizes for m in masks_of_size(n, s)]
    reps = [[]]
    frontier = [[]]
    level = 0
    while frontier:
        new = []
        for fam in frontier:
            for b in big:
                if b in fam or any(popcount(b & c) > 2 for c in fam):
                    continue
                child = sorted(fam + [b])
                if not usage_ok(n, child, caps):
                    continue
                inv = invariant(n, child)
                if any(inv == inv2 and is_isomorphic(n, child, f2) for inv2, f2 in new):
                    continue
                new.append((inv, child))
        frontier = [f for _, f in new]
        level += 1
        if verbose:
            print(f"  phase1 level {level}: {len(frontier)} classes", flush=True)
        reps.extend(frontier)
    return reps


# ----------------------------------------------------------------------------------------------
# lines: maximum / all sets of pairwise <=1-intersecting lines (blocks of F or uncovered triples)
# ----------------------------------------------------------------------------------------------
def uncovered_triples(n, F):
    return [t for t in masks_of_size(n, 3) if not any(b & t == t for b in F)]


def line_search(n, F, caps, min_size=None):
    """min_size None: return (max |L|, one maximum L).  Else: return the list of ALL line sets L with
    |L| >= min_size (each exactly once)."""
    cands = list(F) + uncovered_triples(n, F)
    m = len(cands)
    comp = []
    for i in range(m):
        x = 0
        for j in range(m):
            if j != i and popcount(cands[i] & cands[j]) <= 1:
                x |= 1 << j
        comp.append(x)
    pw = [math.comb(popcount(c), 2) for c in cands]
    is3 = [1 if popcount(c) == 3 else 0 for c in cands]
    best = [0, []]; found = []; chosen = []

    def rec(avail, pairs, n3):
        k = len(chosen)
        if min_size is None:
            if k > best[0]:
                best[0] = k; best[1] = list(chosen)
            if k + popcount(avail) <= best[0]:
                return
        else:
            if k >= min_size:
                found.append(list(chosen))
            if k + popcount(avail) < min_size:
                return
        x = avail
        while x:
            low = x & -x; i = low.bit_length() - 1; x ^= low
            if pairs + pw[i] > caps.capL or n3 + is3[i] > caps.t3n:
                continue
            chosen.append(cands[i])
            rec(x & comp[i], pairs + pw[i], n3 + is3[i])
            chosen.pop()
    rec((1 << m) - 1, 0, 0)
    return (best[0], best[1]) if min_size is None else found


# ----------------------------------------------------------------------------------------------
# phase 2: orderly extension of a fixed family R by 4-blocks under Aut(R)
# ----------------------------------------------------------------------------------------------
def phase2(n, caps, R, aut, need, collect_all=False, quad_filter=None, verbose=False,
           max4_override=None):
    """R: fixed block family (list of bitmasks); aut: list of permutations preserving R (must contain
    the identity).  need = C(n,3) - target: a family F is a *candidate* iff D(F) + ellmax(F) >= need.
    collect_all: return every generated family instead (validation).  quad_filter: optional predicate
    restricting the candidate 4-blocks (must be Aut-invariant)."""
    assert tuple(range(n)) in aut
    quads = [m for m in masks_of_size(n, 4) if all(popcount(m & b) <= 2 for b in R)]
    if quad_filter is not None:
        quads = [m for m in quads if quad_filter(m)]
    M = len(quads)
    idx = {m: i for i, m in enumerate(quads)}
    G = len(aut)
    max4 = caps.max4 if max4_override is None else max4_override
    # permutation table on candidate indices
    T = np.empty((G, M), dtype=np.int16)
    if M:
        pts = np.array([bits(m) for m in quads], dtype=np.int64)
        lookup = np.full(1 << n, -1, dtype=np.int64)
        lookup[np.array(quads)] = np.arange(M)
        for gi, s in enumerate(aut):
            sig = np.array(s, dtype=np.int64)
            codes = (1 << sig[pts]).sum(axis=1)
            T[gi] = lookup[codes]
        assert (T >= 0).all(), "aut does not preserve the candidate set"
    compat = []
    for i in range(M):
        x = 0
        for j in range(M):
            if popcount(quads[i] & quads[j]) <= 2:
                x |= 1 << j
        compat.append(x)
    qpts = [bits(m) for m in quads]
    usage = point_usage(n, R)
    deg4 = [0] * n
    D0 = deficit(R)
    D_min = need - caps.ellmax
    F4 = []
    results = []
    stats = {'nodes': 0, 'evals': 0, 'M': M, 'G': G, 'D0': D0}

    def canonical():
        if G == 1:
            return True
        img = T[:, F4]
        rowmin = img.min(axis=1)
        f0 = F4[0]
        if (rowmin < f0).any():
            return False
        sel = np.nonzero(rowmin == f0)[0]
        if len(sel) <= 1:
            return True
        sub = np.sort(img[sel], axis=1)
        Farr = np.array(F4, dtype=np.int16)
        diff = sub != Farr
        has = diff.any(axis=1)
        if not has.any():
            return True
        first = diff.argmax(axis=1)
        vals = sub[np.arange(len(sel)), first]
        return not (has & (vals < Farr[first])).any()

    def evaluate(D):
        F = list(R) + [quads[i] for i in F4]
        ell, L = line_search(n, F, caps)
        if D + ell >= need:
            results.append({'blocks': [bits(b) for b in F], 'D': D, 'ellmax': ell,
                            'count': math.comb(n, 3) - D - ell, 'example_lines': [bits(l) for l in L]})

    def dfs(last, avail, D):
        stats['nodes'] += 1
        if collect_all:
            results.append([bits(b) for b in list(R) + [quads[i] for i in F4]])
        elif D >= D_min:
            stats['evals'] += 1
            evaluate(D)
        hi = avail & ~((1 << (last + 1)) - 1) if last >= 0 else avail
        if not collect_all:
            cnt = popcount(hi)
            slack = sum(min((caps.capD - usage[p]) // 3, max4 - deg4[p]) for p in range(n))
            add = min(cnt, slack // 4)
            if D + 3 * add + caps.ellmax < need:
                return
        x = hi
        while x:
            low = x & -x; i = low.bit_length() - 1; x ^= low
            ps = qpts[i]
            if any(usage[p] + 3 > caps.capD or deg4[p] >= max4 for p in ps):
                continue
            F4.append(i)
            if canonical():
                for p in ps:
                    usage[p] += 3; deg4[p] += 1
                dfs(i, avail & compat[i], D + 3)
                for p in ps:
                    usage[p] -= 3; deg4[p] -= 1
            F4.pop()

    dfs(-1, (1 << M) - 1, D0)
    return results, stats


def root_bound(n, caps, R, need, quad_filter=None):
    """Upper bound for D of any extension of R by 4-blocks (same bound as the DFS uses at the root)."""
    quads = [m for m in masks_of_size(n, 4) if all(popcount(m & b) <= 2 for b in R)]
    if quad_filter is not None:
        quads = [m for m in quads if quad_filter(m)]
    usage = point_usage(n, R)
    slack = sum(min((caps.capD - usage[p]) // 3, caps.max4) for p in range(n))
    add = min(len(quads), slack // 4)
    return deficit(R) + 3 * add, len(quads)


# ----------------------------------------------------------------------------------------------
# hereditary checks on a linear space (points, lines pairwise <=1-intersecting)
# ----------------------------------------------------------------------------------------------
def hereditary_violations(points, lines, t3=None, check83=True, first_only=True):
    """Return a list of violations (kind, subset, info):
       'SG'  : subset S, not inside one line, all of whose pairs lie on lines with >= 3 points of S
               (violates Sylvester-Gallai for the real point set S; includes the Fano plane);
       '8_3' : 8-subset carrying >= 8 three-point lines (Moebius-Kantor (8_3), not realisable);
       'T3'  : r-subset carrying more than t3(r) three-point lines (orchard numbers, cited)."""
    pts = sorted(points)
    lines = [frozenset(l) for l in lines]
    viol = []
    for r in range(3, len(pts) + 1):
        for S in itertools.combinations(pts, r):
            Ss = frozenset(S)
            if any(Ss <= l for l in lines):
                continue
            sub = [l & Ss for l in lines if len(l & Ss) >= 3]
            covered = sum(len(x) * (len(x) - 1) // 2 for x in sub)
            n3 = sum(1 for x in sub if len(x) == 3)
            if covered == r * (r - 1) // 2:
                viol.append(('SG', S, [sorted(x) for x in sub]))
            if check83 and r == 8 and n3 >= 8:
                viol.append(('8_3', S, n3))
            if t3 is not None and r in t3 and n3 > t3[r]:
                viol.append(('T3', S, n3))
            if viol and first_only:
                return viol
    return viol


def derived_lines(n, F, p):
    return [frozenset(bits(b & ~(1 << p))) for b in F if b >> p & 1]


def check_structure(n, F, t3=None, check83=True):
    """Hereditary checks of the derived structure at every point.  Returns {p: violations}."""
    rep = {}
    for p in range(n):
        v = hereditary_violations([q for q in range(n) if q != p], derived_lines(n, F, p), t3, check83)
        if v:
            rep[p] = v
    return rep


# ----------------------------------------------------------------------------------------------
# Miquel closure (classical theorem of the real inversive plane; optional extra filter)
# ----------------------------------------------------------------------------------------------
CUBE_FACES = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 4, 5), (2, 3, 6, 7), (0, 2, 4, 6), (1, 3, 5, 7)]


def miquel_violations(n, F, first_only=True):
    """Cubes (8 distinct points, 6 faces) such that 5 faces lie in 5 *distinct* blocks of F and the
    6th face lies in no block.  By Miquel's theorem such a structure is not realisable."""
    blk = {}
    for b in F:
        for q in itertools.combinations(bits(b), 4):
            blk[mask(q)] = b
    viol = []
    for S in itertools.combinations(range(n), 8):
        s0, rest = S[0], S[1:]
        for perm in itertools.permutations(rest):
            lab = (s0,) + perm
            faces = [mask(lab[i] for i in f) for f in CUBE_FACES]
            got = [blk.get(f) for f in faces]
            missing = [i for i, g in enumerate(got) if g is None]
            if len(missing) != 1:
                continue
            present = [g for g in got if g is not None]
            if len(set(present)) != 5:
                continue
            viol.append({'labelling': lab, 'faces': [bits(f) for f in faces], 'missing_face': missing[0]})
            if first_only:
                return viol
    return viol


def count_circles(n, F, L):
    return math.comb(n, 3) - deficit(F) - len(L)
