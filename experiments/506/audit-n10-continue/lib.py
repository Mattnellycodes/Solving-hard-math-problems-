"""Audit of n10-continue (Erdős #506, n = 10): own utilities, written from scratch.

Conventions.  Points 0..n-1; a block/line is an int bitmask.  A *structure* on P (n = 10) is
(F, L): F = rich blocks (size >= 4, any two share <= 2 points, so every triple is in <= 1 block),
L = lines (size >= 3; a line of size >= 4 is a member of F; a 3-line is a triple not inside any
block of F; lines pairwise share <= 1 point).  circles(P) = C(n,3) - D(F) - |L|,
D(F) = sum (C(|B|,3) - 1).

Möbius structure M on P u {inf} (inf = index n): blocks of size >= 4 are the circles of F (F - L)
and the lines extended by inf.  Inverting about q in P u {inf} gives a real linear space on the
other n points whose rich lines are B - q for the blocks B of M through q (size >= 4).
"""
from itertools import combinations, permutations
from math import comb

INF = 10


def popcount(x):
    return bin(x).count("1")


def bits(x):
    out = []
    while x:
        b = x & -x
        out.append(b.bit_length() - 1)
        x ^= b
    return out


def mask(pts):
    m = 0
    for p in pts:
        m |= 1 << p
    return m


def apply_perm(blocks, perm):
    """perm: tuple, point i -> perm[i].  Returns sorted tuple of relabelled masks."""
    out = []
    for b in blocks:
        m = 0
        for p in bits(b):
            m |= 1 << perm[p]
        out.append(m)
    return tuple(sorted(out))


# ---------------------------------------------------------------- canonical form (own I/R)

def _refine(blocks, n, cells):
    """Refine an ordered partition (list of lists of points) by block signatures until stable."""
    while True:
        cell_of = {}
        for ci, c in enumerate(cells):
            for p in c:
                cell_of[p] = ci
        sig = {}
        for p in range(n):
            s = []
            for b in blocks:
                if b >> p & 1:
                    s.append((popcount(b), tuple(sorted(cell_of[q] for q in bits(b) if q != p))))
            sig[p] = tuple(sorted(s))
        new = []
        for c in cells:
            groups = {}
            for p in c:
                groups.setdefault(sig[p], []).append(p)
            for key in sorted(groups):
                new.append(groups[key])
        if len(new) == len(cells):
            return new
        cells = new


def canon(blocks, n):
    """Canonical form: (canonical sorted tuple of masks, canonical perm, |Aut|).
    Individualisation-refinement search over all leaves, with twin pruning: points with identical
    incidence sets ("twins") are interchanged by an automorphism, so only one twin per cell is
    individualised.  |Aut| = (#leaves attaining the minimum) * prod (twin-class sizes)!."""
    blocks = list(blocks)
    best = None
    best_perm = None
    count = 0
    inc = [frozenset(b for b in blocks if b >> p & 1) for p in range(n)]
    twins = {}
    for p in range(n):
        twins.setdefault(inc[p], []).append(p)
    tw = 1
    for t in twins.values():
        for i in range(2, len(t) + 1):
            tw *= i

    def rec(cells):
        nonlocal best, best_perm, count
        cells = _refine(blocks, n, cells)
        for ci, c in enumerate(cells):
            if len(c) > 1:
                seen = set()
                for p in c:
                    if inc[p] in seen:
                        continue
                    seen.add(inc[p])
                    rest = [q for q in c if q != p]
                    rec(cells[:ci] + [[p], rest] + cells[ci + 1:])
                return
        perm = [0] * n
        for i, c in enumerate(cells):
            perm[c[0]] = i
        img = apply_perm(blocks, perm)
        if best is None or img < best:
            best, best_perm, count = img, tuple(perm), 1
        elif img == best:
            count += 1

    rec([list(range(n))])
    return best, best_perm, count * tw


def automorphisms(blocks, n):
    """All permutations of 0..n-1 mapping the block set onto itself (backtracking)."""
    bset = set(blocks)
    blocks = list(blocks)
    inc = [[b for b in blocks if b >> p & 1] for p in range(n)]
    # pair-count and degree invariants for pruning
    deg = [tuple(sorted(popcount(b) for b in inc[p])) for p in range(n)]
    pair = {}
    for p in range(n):
        for q in range(n):
            pair[p, q] = tuple(sorted(popcount(b) for b in inc[p] if b >> q & 1))
    out = []
    img = [-1] * n
    used = [False] * n

    def rec(i):
        if i == n:
            if apply_perm(blocks, img) == tuple(sorted(bset)):
                out.append(tuple(img))
            return
        for t in range(n):
            if used[t] or deg[t] != deg[i]:
                continue
            ok = True
            for j in range(i):
                if pair[i, j] != pair[t, img[j]]:
                    ok = False
                    break
            if not ok:
                continue
            img[i] = t
            used[t] = True
            rec(i + 1)
            used[t] = False
            img[i] = -1

    rec(0)
    return out


# ---------------------------------------------------------------- structures

def D_of(F):
    return sum(comb(popcount(b), 3) - 1 for b in F)


def circles_count(F, L, n=10):
    return comb(n, 3) - D_of(F) - len(L)


def check_structure(F, L, n=10):
    """Basic combinatorial validity of (F, L)."""
    F = list(F)
    for a, b in combinations(F, 2):
        assert popcount(a & b) <= 2, "two rich blocks share 3 points"
    for b in F:
        assert popcount(b) >= 4
    for l in L:
        k = popcount(l)
        assert k >= 3
        if k >= 4:
            assert l in F, "rich line not a block"
        else:
            assert all(popcount(l & b) < 3 for b in F), "3-line inside a rich block"
    for a, b in combinations(L, 2):
        assert popcount(a & b) <= 1, "two lines share 2 points"
    return True


def mobius_blocks(F, L, n=10):
    """Blocks of size >= 4 of the Möbius structure on P u {inf}."""
    Ls = set(L)
    M = [b for b in F if b not in Ls]
    M += [l | (1 << n) for l in L]
    return M


def derived(M, q):
    """Rich lines (masks on the remaining points, point indices unchanged) of the real linear
    space obtained by inverting the Möbius structure M about q."""
    return [b & ~(1 << q) for b in M if b >> q & 1]


def all_derived(F, L, n=10):
    """List of (q, points_mask, lines) for q in P u {inf}."""
    M = mobius_blocks(F, L, n)
    full = (1 << (n + 1)) - 1
    return [(q, full & ~(1 << q), derived(M, q)) for q in range(n + 1)]


# ---------------------------------------------------------------- real-linear-space filters

def sg_pairs(lines):
    return sum(comb(popcount(l), 2) for l in lines)


def sg_cap_ok(lines, npts):
    """Sylvester-Gallai with o(m) >= 1: rich lines cover <= C(m,2) - 1 pairs."""
    return sg_pairs(lines) <= comb(npts, 2) - 1


def submasks(pm):
    """All submasks of pm (including 0)."""
    s = pm
    while True:
        yield s
        if s == 0:
            return
        s = (s - 1) & pm


def hereditary_sg_violation(lines, pm):
    """Return a subset S (mask) of pm, |S| >= 3, not inside one line, all of whose pairs are
    covered by lines restricted to S; or None."""
    for S in submasks(pm):
        k = popcount(S)
        if k < 3:
            continue
        covered = 0
        inside = False
        for l in lines:
            r = l & S
            c = popcount(r)
            if c >= 3:
                covered += comb(c, 2)
            if r == S:
                inside = True
        if not inside and covered == comb(k, 2):
            return S
    return None


T3CAP_THEOREM = {3: 1, 4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10, 10: 13}
# 3..6: combinatorial packing bounds; 7: Fano is not real (SG); 8: (8_3) is not real;
# 9: every 11-packing contains an (8_3) (checked in sts9_check.py); 10: each point on <= 4 lines.
T3CAP_CITED = {**T3CAP_THEOREM, 9: 10, 10: 12}


def orchard_violation(lines, pm, caps=T3CAP_THEOREM):
    """Return (S, count) if some subset S of pm has more lines with exactly 3 points of S than
    caps[|S|] allows; else None."""
    for S in submasks(pm):
        k = popcount(S)
        if k < 7:
            continue
        c = sum(1 for l in lines if popcount(l & S) == 3)
        if c > caps[k]:
            return (S, c)
    return None


def mk_violation(lines, pm):
    """An 8-subset with >= 8 three-point lines (=> contains a Möbius-Kantor (8_3))."""
    pts = bits(pm)
    for S8 in combinations(pts, 8):
        S = mask(S8)
        c = sum(1 for l in lines if popcount(l & S) == 3)
        if c >= 8:
            return S
    return None


# ---------------------------------------------------------------- Miquel

def miquel_violation(M, n=10):
    """M: blocks (size >= 4) of a Möbius structure on n+1 points.  Return a cube (8 points
    labelled by the faces) with exactly 5 co-blocked faces, or None.  Faces of the cube: the two
    x-faces Q0, Q1 (disjoint quadruples) are co-blocked; the four side faces are determined by a
    bijection Q0 -> Q1 and a 4-cycle on Q0."""
    quads = set()
    for b in M:
        for q in combinations(bits(b), 4):
            quads.add(mask(q))
    quads = sorted(quads)
    qset = set(quads)

    def coblocked(a, b, c, d):
        return mask((a, b, c, d)) in qset

    cycles = [(0, 1, 2, 3), (0, 1, 3, 2), (0, 2, 1, 3)]
    for i, Q0 in enumerate(quads):
        p = bits(Q0)
        for Q1 in quads[i + 1:]:
            if Q0 & Q1:
                continue
            q = bits(Q1)
            for cyc in cycles:
                a, b, c, d = (p[t] for t in cyc)
                for perm in permutations(range(4)):
                    a1, b1, c1, d1 = (q[t] for t in perm)
                    cnt = 2
                    cnt += coblocked(a, b, a1, b1)
                    cnt += coblocked(c, d, c1, d1)
                    cnt += coblocked(a, d, a1, d1)
                    cnt += coblocked(b, c, b1, c1)
                    if cnt == 5:
                        return ((a, b, c, d), (a1, b1, c1, d1))
    return None


# ---------------------------------------------------------------- all filters

def run_filters(F, L, n=10, caps=T3CAP_THEOREM, with_mk=True):
    """Return list of (filter name, q, witness) violations for the structure (F, L)."""
    viol = []
    for q, pm, lines in all_derived(F, L, n):
        if not sg_cap_ok(lines, n):
            viol.append(("sg_cap", q, sg_pairs(lines)))
        w = hereditary_sg_violation(lines, pm)
        if w is not None:
            viol.append(("hereditary_sg", q, w))
        if with_mk:
            w = mk_violation(lines, pm)
            if w is not None:
                viol.append(("mk_8_3", q, w))
        w = orchard_violation(lines, pm, caps)
        if w is not None:
            viol.append(("orchard", q, w))
    w = miquel_violation(mobius_blocks(F, L, n), n)
    if w is not None:
        viol.append(("miquel", None, w))
    return viol


# ---------------------------------------------------------------- exact geometry -> structure

def structure_from_points(points):
    """Exact (Fraction) coordinates -> (F, L) masks.  Own implementation."""
    from fractions import Fraction as Fr
    pts = [(Fr(x), Fr(y)) for x, y in points]
    n = len(pts)
    assert len(set(pts)) == n
    groups = {}
    for i, j, k in combinations(range(n), 3):
        (ax, ay), (bx, by), (cx, cy) = pts[i], pts[j], pts[k]
        d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
        if d == 0:
            # line through i, j: normalised (a, b, c) with a x + b y = c
            a, b = by - ay, ax - bx
            c = a * ax + b * ay
            if a != 0:
                b, c, a = b / a, c / a, Fr(1)
            else:
                c, b = c / b, Fr(1)
            key = ("L", a, b, c)
        else:
            a2, b2, c2 = ax * ax + ay * ay, bx * bx + by * by, cx * cx + cy * cy
            ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d
            uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d
            key = ("C", ux, uy, (ax - ux) ** 2 + (ay - uy) ** 2)
        groups.setdefault(key, 0)
        groups[key] |= mask((i, j, k))
    F = sorted(m for m in groups.values() if popcount(m) >= 4)
    L = sorted(m for key, m in groups.items() if key[0] == "L")
    return F, L


if __name__ == "__main__":
    # self-tests
    # Fano plane: hereditary SG violation; (8_3): mk violation
    fano = [mask(t) for t in [(0, 1, 2), (0, 3, 4), (0, 5, 6), (1, 3, 5), (1, 4, 6), (2, 3, 6), (2, 4, 5)]]
    assert hereditary_sg_violation(fano, 127) == 127
    mk = [mask(((i) % 8, (i + 1) % 8, (i + 3) % 8)) for i in range(8)]
    assert mk_violation(mk, 255) == 255
    assert hereditary_sg_violation(mk, 255) is None
    # antipodal configuration: 9 points on the unit circle (4 antipodal pairs + 1 more)
    # and the centre.  Count must be 33 and all filters pass.
    from fractions import Fraction as Fr
    circ = [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(3, 5), Fr(4, 5)), (Fr(-4, 5), Fr(3, 5))]
    pts = []
    for x, y in circ:
        pts += [(x, y), (-x, -y)]
    pts += [(Fr(0), Fr(0)), (Fr(5, 13), Fr(12, 13))]   # centre + a 9th point on the circle
    F, L = structure_from_points(pts)
    check_structure(F, L)
    print("antipodal: |F| =", len(F), "sizes", sorted(popcount(b) for b in F), "|L| =", len(L),
          "circles =", circles_count(F, L))
    assert circles_count(F, L) == 33
    v = run_filters(F, L)
    print("antipodal filter violations:", v)
    assert not v
    # Miquel: the cube (8 points of a real configuration) is closed; removing one face is not
    cube = [mask(q) for q in [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 4, 5), (2, 3, 6, 7), (0, 3, 4, 7), (1, 2, 5, 6)]]
    assert miquel_violation(cube, 7) is None
    assert miquel_violation(cube[:5], 7) is not None
    # canonical form / automorphisms sanity
    cf, perm, na = canon(fano, 7)
    assert na == 168, na
    assert len(automorphisms(fano, 7)) == 168
    cf, perm, na = canon(mk, 8)
    assert na == 48 and len(automorphisms(mk, 8)) == 48
    import random
    random.seed(1)
    for trial in range(40):
        nn = random.choice([6, 7, 8, 9, 10])
        fam = []
        for _ in range(random.randint(1, 6)):
            k = random.choice([3, 3, 4, 4, 5, 6])
            b = mask(random.sample(range(nn), k))
            if all(popcount(b & c) <= 2 for c in fam):
                fam.append(b)
        na = canon(fam, nn)[2]
        nb = len(automorphisms(fam, nn))
        assert na == nb, (fam, na, nb)
        perm = list(range(nn)); random.shuffle(perm)
        assert canon(apply_perm(fam, perm), nn)[0] == canon(fam, nn)[0]
    assert canon([mask(range(9))], 10)[2] == 362880
    print("self-tests passed")
