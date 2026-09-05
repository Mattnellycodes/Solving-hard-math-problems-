"""Hostile-referee checks of the auxiliary lemmas quoted in REPORT.md (not all are load-bearing
for c(6); see the audit report for which ones are).

 1. (8_3): every family of 8 triples on 8 points pairwise sharing <= 1 point is isomorphic to the
    Mobius-Kantor configuration MK = {i, i+1, i+3 mod 8}  (my own brute force, my own canonical form),
    and MK has no real projective realisation -- proved here with a DIFFERENT frame {0,2,4,6} from the
    theory agent's {0,1,2,5}:  p0=(1,0,0), p2=(0,1,0), p4=(0,0,1), p6=(1,1,1); the line structure forces
    p1=(0,1,v), p3=(1,1,w), p5=(u,1,1), p7=(1,r,0) (alternatives coincide with frame points), and the
    remaining four incidences give a system whose real solutions we determine exactly.
 2. Fano plane not realisable over R: all 21 pairs covered by 7 three-point lines -> no ordinary line,
    contradicting Sylvester-Gallai.  (Combinatorial check that the 7-line 3-uniform pairwise <=1
    systems on 7 points covering all pairs are all Fano.)
 3. Hereditary Sylvester-Gallai applied to the unique n = 6 candidate: derived structures at every point
    and the line structure at infinity.  (We expect NO violation -- i.e. hereditary SG does not kill the
    candidate; only the Angle Lemma does.)
 4. Lemma A (block of size n-1) recomputed by brute force for n = 5..8: minimum combinatorial count of a
    structure containing an (n-1)-block, with the block a circle or a line.
"""
import itertools
from math import comb
import sympy as sp


def canon_triples(F, n):
    return min(tuple(sorted(tuple(sorted(s[i] for i in t)) for t in F)) for s in itertools.permutations(range(n)))


def all_partial_linear_spaces(n, k, nlines):
    """all families of nlines k-subsets of [n] pairwise sharing <= 1 point, up to isomorphism."""
    triples = [frozenset(t) for t in itertools.combinations(range(n), k)]
    classes = set()

    def rec(F, start):
        if len(F) == nlines:
            classes.add(canon_triples(F, n)); return
        for j in range(start, len(triples)):
            t = triples[j]
            if all(len(t & u) <= 1 for u in F):
                rec(F + [t], j + 1)
    rec([], 0)
    return classes


def check_83():
    classes = all_partial_linear_spaces(8, 3, 8)
    MK = [frozenset({i, (i + 1) % 8, (i + 3) % 8}) for i in range(8)]
    print(f"[8_3] isomorphism classes of 8 triples on 8 points pairwise sharing <=1 point: {len(classes)}; "
          f"MK among them: {canon_triples(MK, 8) in classes}")
    v, w, u, r = sp.symbols('v w u r')
    p = {0: sp.Matrix([1, 0, 0]), 2: sp.Matrix([0, 1, 0]), 4: sp.Matrix([0, 0, 1]), 6: sp.Matrix([1, 1, 1]),
         1: sp.Matrix([0, 1, v]), 3: sp.Matrix([1, 1, w]), 5: sp.Matrix([u, 1, 1]), 7: sp.Matrix([1, r, 0])}
    # sanity: the frame {0,2,4,6} has no three points on an MK line
    assert not any(len(L & {0, 2, 4, 6}) >= 3 for L in MK)
    # sanity: p1 on line(2,4)={x=0}, p3 on line(4,6)={x=y}, p5 on line(6,0)={y=z}, p7 on line(0,2)={z=0}
    eqs = [sp.expand(sp.Matrix.hstack(*[p[i] for i in sorted(L)]).det()) for L in MK]
    eqs = [e for e in eqs if e != 0]
    print("[8_3] remaining incidence equations:", eqs)
    G = sp.groebner(eqs, v, w, u, r, order='lex')
    print("[8_3] lex Groebner basis:", [sp.factor(g) for g in G])
    uni = [g for g in G if len(g.free_symbols) == 1]
    for g in uni:
        x = list(g.free_symbols)[0]
        print(f"[8_3] univariate {x}: {sp.factor(g)}  real roots: {sp.real_roots(sp.Poly(g, x))}")
    # distinctness also needed: v != 0 (p1 != p2), w != 1 (p3 != p6), u != 1 (p5 != p6), r != 0 (p7 != p0)
    sols = sp.solve(eqs, [v, w, u, r], dict=True)
    print("[8_3] all complex solutions:", sols)
    real_sols = [s for s in sols if all(sp.im(sp.nsimplify(val)).equals(0) for val in s.values())]
    print("[8_3] real solutions:", real_sols)


def check_fano():
    classes = all_partial_linear_spaces(7, 3, 7)
    full = [c for c in classes if len({frozenset(pq) for t in c for pq in itertools.combinations(t, 2)}) == 21]
    print(f"[Fano] classes of 7 triples on 7 points pairwise sharing <=1 point: {len(classes)}; "
          f"classes covering all 21 pairs (Fano): {len(full)} -> such a derived structure has no ordinary "
          f"line, contradicting Sylvester-Gallai")


def sg_closed_subsets(points, lines):
    """subsets S (|S|>=3) not inside one line such that every pair of S is in a line meeting S in >=3 pts."""
    viol = []
    pts = sorted(points)
    for k in range(3, len(pts) + 1):
        for S in itertools.combinations(pts, k):
            Ss = set(S)
            if any(Ss <= l for l in lines):
                continue
            cov = set()
            for l in lines:
                m = l & Ss
                if len(m) >= 3:
                    cov |= {frozenset(pq) for pq in itertools.combinations(m, 2)}
            if len(cov) == comb(k, 2):
                viol.append(S)
    return viol


def check_hereditary_n6():
    blocks = [frozenset(b) for b in ([0, 1, 2, 3], [0, 1, 4, 5], [2, 3, 4, 5])]
    lines = [frozenset(l) for l in ([0, 2, 4], [0, 3, 5], [1, 2, 5], [1, 3, 4])]
    for pnt in range(6):
        derived = [b - {pnt} for b in blocks if pnt in b]
        v = sg_closed_subsets(set(range(6)) - {pnt}, derived)
        print(f"[hereditary SG, n=6 candidate] point {pnt}: derived lines {[sorted(x) for x in derived]} "
              f"violations: {v}")
    v = sg_closed_subsets(set(range(6)), lines)
    print(f"[hereditary SG, n=6 candidate] lines at infinity {[sorted(x) for x in lines]} violations: {v}")


def check_lemma_A():
    for n in range(5, 9):
        pts = range(n)
        B = frozenset(range(n - 1)); x = n - 1
        N3 = comb(n, 3)
        # the only blocks: B and {x,a,b}; lines: either B (line) alone, or 3-lines {x,a,b} with disjoint pairs
        # brute force over line sets for both cases
        triples = [frozenset((x, a, b)) for a, b in itertools.combinations(range(n - 1), 2)]
        best_circle = None; best_line = None
        # case B circle: lines subset of triples pairwise sharing <=1
        def rec(L, start):
            nonlocal best_circle
            best_circle = max(best_circle or 0, len(L))
            for j in range(start, len(triples)):
                if all(len(triples[j] & T) <= 1 for T in L):
                    rec(L + [triples[j]], j + 1)
        rec([], 0)
        D = comb(n - 1, 3) - 1
        circ_case = N3 - D - best_circle
        line_case = N3 - D - 1
        f = comb(n - 1, 2) + 1 - (n - 1) // 2
        print(f"[Lemma A] n={n}: (n-1)-block a circle -> min circles {circ_case} (ell_max={best_circle}); "
              f"a line -> {line_case}; f(n)={f}; Lemma A claims min = f(n): {min(circ_case, line_case) == f}")


if __name__ == "__main__":
    check_83()
    check_fano()
    check_hereditary_n6()
    check_lemma_A()
