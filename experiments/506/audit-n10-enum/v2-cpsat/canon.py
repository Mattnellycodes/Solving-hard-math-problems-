"""Canonical labelling of set systems (hypergraphs) on n points, written from scratch.

canon(n, blocks) -> (code, relabelled_blocks, aut_order, automorphisms)
  code: tuple of sorted bitmasks of the relabelled blocks (isomorphism invariant, complete)
Method: individualisation-refinement search tree; canonical code = min over leaves; pruning of
children by orbits of automorphisms found so far (sound: an automorphism fixing the individualised
prefix maps subtrees to subtrees with identical leaf codes).
"""
import itertools, random

def refine(n, blocks, cells):
    # cells: list of lists (ordered partition). returns refined ordered partition (equivariant).
    cell_of = [0]*n
    while True:
        for i, c in enumerate(cells):
            for v in c: cell_of[v] = i
        sig = {}
        for v in range(n):
            s = []
            for B in blocks:
                if v in B:
                    s.append((len(B), tuple(sorted(cell_of[u] for u in B if u != v))))
            sig[v] = tuple(sorted(s))
        new = []
        for c in cells:
            if len(c) == 1:
                new.append(c); continue
            groups = {}
            for v in c:
                groups.setdefault(sig[v], []).append(v)
            for k in sorted(groups):
                new.append(groups[k])
        if len(new) == len(cells):
            return new
        cells = new

def canon(n, blocks, want_auts=False, init_cells=None):
    blocks = [frozenset(B) for B in blocks]
    best = [None, None]   # code, labelling
    first = [None, None]
    auts = []
    def leaf(cells):
        lab = [0]*n
        for i, c in enumerate(cells): lab[c[0]] = i
        code = tuple(sorted(sum(1 << lab[v] for v in B) for B in blocks))
        for ref in (first, best):
            if ref[0] is not None and ref[0] == code and ref[1] != lab:
                inv = [0]*n
                for v in range(n): inv[ref[1][v]] = v
                sigma = tuple(inv[lab[v]] for v in range(n))   # sigma maps v -> point with same label in ref
                auts.append(sigma)
        if first[0] is None:
            first[0], first[1] = code, lab
        if best[0] is None or code < best[0]:
            best[0], best[1] = code, list(lab)
    def search(cells, fixed):
        cells = refine(n, blocks, cells)
        if len(cells) == n:
            leaf(cells); return
        # first non-singleton cell
        idx = next(i for i, c in enumerate(cells) if len(c) > 1)
        cell = cells[idx]
        done = []
        for v in cell:
            # prune: is v in the orbit of an explored child under automorphisms fixing 'fixed' pointwise?
            gens = [a for a in auts if all(a[f] == f for f in fixed)]
            if gens and done:
                orb = set(done); frontier = list(done)
                while frontier:
                    u = frontier.pop()
                    for a in gens:
                        w = a[u]
                        if w not in orb:
                            orb.add(w); frontier.append(w)
                if v in orb:
                    continue
            newcells = cells[:idx] + [[v], [u for u in cell if u != v]] + cells[idx+1:]
            search(newcells, fixed + [v])
            done.append(v)
    search([list(c) for c in init_cells] if init_cells else [list(range(n))], [])
    lab = best[1]
    rel = sorted(tuple(sorted(lab[v] for v in B)) for B in blocks)
    return best[0], rel, lab, auts

def aut_group_order(n, blocks):
    """order of Aut by closing the generators found (brute-force closure; fine for small groups)"""
    _, _, _, auts = canon(n, blocks)
    gens = set(auts)
    ident = tuple(range(n))
    group = {ident}
    frontier = [ident]
    while frontier:
        g = frontier.pop()
        for a in gens:
            h = tuple(a[g[i]] for i in range(n))
            if h not in group:
                group.add(h); frontier.append(h)
    return len(group), group

def brute_canon(n, blocks):
    blocks = [frozenset(B) for B in blocks]
    best = None; cnt = 0
    for p in itertools.permutations(range(n)):
        code = tuple(sorted(sum(1 << p[v] for v in B) for B in blocks))
        if best is None or code < best:
            best = code; cnt = 1
        elif code == best:
            cnt += 1
    return best, cnt   # cnt = |Aut|

if __name__ == '__main__':
    random.seed(1)
    # validation 1: the equivalence relation induced by canon() coincides with true isomorphism
    # (brute force over S_n) on a pool of random systems, n <= 7; |Aut| agrees too.
    pool = []
    for trial in range(300):
        n = random.choice([6, 7])
        k = random.randint(1, 7)
        blocks = list({frozenset(random.sample(range(n), random.choice([3, 4, 4, 5]))) for _ in range(k)})
        p = list(range(n)); random.shuffle(p)
        blocks2 = [frozenset(p[v] for v in B) for B in blocks]   # an isomorphic copy
        pool.append((n, blocks)); pool.append((n, blocks2))
    mine = {}; brute = {}; bad = 0
    for n, blocks in pool:
        c1, rel, lab, auts = canon(n, blocks)
        c2, autn = brute_canon(n, blocks)
        o, _ = aut_group_order(n, blocks)
        if o != autn: bad += 1; print('AUT MISMATCH', n, blocks, o, autn)
        mine.setdefault((n, c1), set()).add((n, c2))
        brute.setdefault((n, c2), set()).add((n, c1))
    ok = all(len(v) == 1 for v in mine.values()) and all(len(v) == 1 for v in brute.values())
    print('validation 1: pool', len(pool), 'my classes', len(mine), 'brute classes', len(brute),
          'bijection:', ok, 'aut mismatches:', bad)
    # validation 2: invariance under random relabelling, n = 10, structured examples
    bad = 0
    for trial in range(200):
        n = 10
        blocks = set()
        while len(blocks) < random.randint(5, 22):
            s = random.choice([4, 4, 4, 5])
            B = frozenset(random.sample(range(n), s))
            if all(len(B & C) <= 2 for C in blocks): blocks.add(B)
        blocks = list(blocks)
        c1 = canon(n, blocks)[0]
        p = list(range(n)); random.shuffle(p)
        c2 = canon(n, [frozenset(p[v] for v in B) for B in blocks])[0]
        if c1 != c2: bad += 1; print('NOT INVARIANT', blocks)
    print('validation 2 (relabelling invariance, n=10): failures =', bad)
    import time
    t = time.time()
    print('Aut of single 5-block on 10 points:', aut_group_order(10, [range(5)])[0], '(expect 14400)')
    print('Aut of two disjoint 5-blocks:', aut_group_order(10, [range(5), range(5,10)])[0], '(expect 28800)')
    print('Aut of 6-block:', aut_group_order(10, [range(6)])[0], '(expect 17280)')
    print('time', time.time() - t)
