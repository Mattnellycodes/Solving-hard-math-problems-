"""Smith normal form with transforms (auditor's own): for an integer matrix M (m x n) return
(D, V) with D = U M V diagonal (U, V unimodular; U not returned), d_1 | d_2 | ... .
Used to describe {phi in R^n : M phi in Z^m} modulo Z^n:  psi = V^{-1} phi has d_i psi_i in Z for
i <= rank, psi_i free for i > rank; and phi == phi' (mod Z^n) iff psi == psi' (mod Z^n)."""
from fractions import Fraction
import random


def snf(M):
    A = [list(map(int, row)) for row in M]
    m, n = len(A), len(A[0])
    V = [[int(i == j) for j in range(n)] for i in range(n)]

    def swap_rows(i, j):
        A[i], A[j] = A[j], A[i]

    def swap_cols(i, j):
        for r in A:
            r[i], r[j] = r[j], r[i]
        for r in V:
            r[i], r[j] = r[j], r[i]

    def add_row(i, j, q):  # row_i -= q row_j
        A[i] = [a - q * b for a, b in zip(A[i], A[j])]

    def add_col(i, j, q):  # col_i -= q col_j
        for r in A:
            r[i] -= q * r[j]
        for r in V:
            r[i] -= q * r[j]

    t = 0
    while t < min(m, n):
        # find smallest nonzero |entry| in submatrix
        best = None
        for i in range(t, m):
            for j in range(t, n):
                if A[i][j] != 0 and (best is None or abs(A[i][j]) < best[0]):
                    best = (abs(A[i][j]), i, j)
        if best is None:
            break
        _, i, j = best
        swap_rows(t, i)
        swap_cols(t, j)
        while True:
            changed = False
            for i in range(t + 1, m):
                if A[i][t] != 0:
                    q = A[i][t] // A[t][t]
                    add_row(i, t, q)
                    if A[i][t] != 0:
                        swap_rows(t, i)
                        changed = True
            for j in range(t + 1, n):
                if A[t][j] != 0:
                    q = A[t][j] // A[t][t]
                    add_col(j, t, q)
                    if A[t][j] != 0:
                        swap_cols(t, j)
                        changed = True
            if changed:
                continue
            # divisibility
            bad = None
            for i in range(t + 1, m):
                for j in range(t + 1, n):
                    if A[i][j] % A[t][t] != 0:
                        bad = i
                        break
                if bad is not None:
                    break
            if bad is None:
                break
            add_row(t, bad, -1)   # row_t += row_bad
        if A[t][t] < 0:
            for r in A:
                r[t] = -r[t]
            for r in V:
                r[t] = -r[t]
        t += 1
    d = [A[i][i] for i in range(min(m, n))]
    return d, V


def det(M):
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    d = Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if A[r][i] != 0), None)
        if p is None:
            return Fraction(0)
        if p != i:
            A[i], A[p] = A[p], A[i]
            d = -d
        d *= A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            A[r] = [a - f * b for a, b in zip(A[r], A[i])]
    return d


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


if __name__ == '__main__':
    random.seed(0)
    for trial in range(300):
        m, n = random.randint(1, 6), random.randint(1, 6)
        M = [[random.randint(-3, 3) for _ in range(n)] for _ in range(m)]
        d, V = snf(M)
        assert abs(det(V)) == 1
        MV = matmul(M, V)
        # MV should have the same column lattice as D up to row ops: check via SNF invariants:
        # verify that MV = U^{-1} D, i.e. columns j > rank of MV are zero and the column lattice
        # of MV[:, :r] is generated with elementary divisors d: simply recompute snf of MV and compare
        d2, _ = snf(MV)
        assert d2 == d, (M, d, d2)
        r = sum(1 for x in d if x != 0)
        for j in range(r, n):
            assert all(row[j] == 0 for row in MV)
        # brute force check: for r == n every solution of M phi in Z^m has denominators dividing
        # e = lcm(d); on the grid (1/e)Z^n / Z^n the solution set must equal V.{psi : d_i psi_i in Z}
        if r == n and n <= 3:
            import itertools, math
            e = 1
            for x in d[:r]:
                e = e * x // math.gcd(e, x)
            if e <= 12:
                sols = set()
                for ks in itertools.product(range(e), repeat=n):
                    phi = [Fraction(k, e) for k in ks]
                    if all(sum(M[i][j] * phi[j] for j in range(n)).denominator == 1 for i in range(m)):
                        sols.add(tuple(phi))
                pred = set()
                for ks in itertools.product(*[range(x) for x in d[:r]]):
                    psi = [Fraction(k, x) for k, x in zip(ks, d[:r])]
                    phi = [sum(V[i][j] * psi[j] for j in range(n)) % 1 for i in range(n)]
                    pred.add(tuple(phi))
                assert pred == sols, (M, d, V, pred, sols)
    print('snf self-test passed')
