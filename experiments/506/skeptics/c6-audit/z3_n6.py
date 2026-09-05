"""Second independent check of the n = 6 combinatorial relaxation, with a different solver (z3, SAT
enumeration via blocking clauses) and a different encoding from the theory agent's CP-SAT model.

Variables: x_B (B rich block, 4 <= |B| <= n-1), y_S (S a line, |S| >= 3).
Constraints: pairwise |B ∩ B'| >= 3 -> not both; |S ∩ S'| >= 2 -> not both; a line of size >= 4 must be
a rich block; a 3-line must not be inside a rich block (a triple inside a rich block is not a block);
optional Sylvester-Gallai pair caps.  We enumerate every labelled solution with
      D + ell = sum (C(|B|,3)-1) x_B + sum y_S  >=  C(n,3) - target
and canonicalise under S_n.  We also compute the exact maximum of D + ell by bisection.
"""
import sys
from itertools import combinations, permutations
from math import comb
import z3

O_STRONG = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6}


def o_val(regime, m):
    return 0 if regime == "none" else (1 if regime == "weak" else O_STRONG[m])


def build(n, regime):
    pts = range(n)
    rich = [frozenset(c) for k in range(4, n) for c in combinations(pts, k)]
    lines = [frozenset(c) for k in range(3, n) for c in combinations(pts, k)]
    x = {B: z3.Bool("x_" + "".join(map(str, sorted(B)))) for B in rich}
    y = {S: z3.Bool("y_" + "".join(map(str, sorted(S)))) for S in lines}
    s = z3.Solver()
    for A, B in combinations(rich, 2):
        if len(A & B) >= 3:
            s.add(z3.Not(z3.And(x[A], x[B])))
    for S, T in combinations(lines, 2):
        if len(S & T) >= 2:
            s.add(z3.Not(z3.And(y[S], y[T])))
    for S in lines:
        if len(S) >= 4:
            s.add(z3.Implies(y[S], x[S]))
        else:
            for B in rich:
                if S <= B:
                    s.add(z3.Not(z3.And(y[S], x[B])))
    cap_d = comb(n - 1, 2) - o_val(regime, n - 1)
    cap_l = comb(n, 2) - o_val(regime, n)
    for p in pts:
        s.add(z3.Sum([z3.If(x[B], comb(len(B) - 1, 2), 0) for B in rich if p in B]) <= cap_d)
    s.add(z3.Sum([z3.If(y[S], comb(len(S), 2), 0) for S in lines]) <= cap_l)
    obj = z3.Sum([z3.If(x[B], comb(len(B), 3) - 1, 0) for B in rich] + [z3.If(y[S], 1, 0) for S in lines])
    return s, x, y, obj


def canon(F, L, n):
    return min((tuple(sorted(tuple(sorted(s[i] for i in B)) for B in F)),
                tuple(sorted(tuple(sorted(s[i] for i in S)) for S in L))) for s in permutations(range(n)))


def main():
    n = int(sys.argv[1]); target = int(sys.argv[2])
    N3 = comb(n, 3)
    for regime in ("none", "weak", "strong"):
        s, x, y, obj = build(n, regime)
        # exact maximum of D + ell
        best = None
        s.push()
        lo, hi = 0, N3 + comb(n, 2)
        while lo < hi:
            mid = (lo + hi + 1) // 2
            s.push(); s.add(obj >= mid)
            if s.check() == z3.sat:
                lo = mid
            else:
                hi = mid - 1
            s.pop()
        best = lo
        s.pop()
        print(f"n={n} regime={regime}: max D+ell = {best}  -> combinatorial minimum circles = {N3 - best}")
        # enumerate all labelled solutions with count <= target
        s.push(); s.add(obj >= N3 - target)
        sols = []
        while s.check() == z3.sat:
            m = s.model()
            F = [B for B in x if z3.is_true(m.eval(x[B], model_completion=True))]
            L = [S for S in y if z3.is_true(m.eval(y[S], model_completion=True))]
            sols.append((F, L))
            lits = [x[B] for B in F] + [y[S] for S in L]
            s.add(z3.Or([z3.Not(v) for v in lits] + [v for v in list(x.values()) + list(y.values()) if v not in lits]))
            if len(sols) > 5000:
                print("   too many solutions, stopping"); break
        s.pop()
        classes = {}
        for F, L in sols:
            classes.setdefault(canon(F, L, n), 0)
            classes[canon(F, L, n)] += 1
        print(f"   labelled solutions with circles <= {target}: {len(sols)}; isomorphism classes: {len(classes)}")
        for (F, L), cnt in sorted(classes.items()):
            D = sum(comb(len(B), 3) - 1 for B in F)
            print(f"   class: blocks={[list(b) for b in F]} lines={[list(s_) for s_ in L]} D={D} ell={len(L)} "
                  f"circles={N3 - D - len(L)} labelled copies={cnt}")


if __name__ == "__main__":
    main()
