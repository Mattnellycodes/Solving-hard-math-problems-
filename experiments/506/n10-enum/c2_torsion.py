#!/usr/bin/env python3
"""c2_torsion.py -- all solutions of the C-2 congruence system  a_i + a_j == b_k + b_l (mod 1)  (twenty equations,
one per 4-block) with all five a's distinct and all five b's distinct.  By the Smith normal form (analyze_C2.py)
the solution group is (R/Z) x Z/10 (rotation x torsion of order 10), so every solution is a rotation of a vector
with coordinates in (1/10) Z.  We compute the kernel of the integer matrix over Z/10 = Z/2 x Z/5 (CRT)."""
import sys, json, itertools
import numpy as np


def system_matrix(rec):
    blocks = [tuple(b) for b in rec['blocks']]
    A, B = [b for b in blocks if len(b) == 5]
    fours = [b for b in blocks if len(b) == 4]
    idx = {p: i for i, p in enumerate(list(A) + list(B))}
    M = np.zeros((len(fours), 10), dtype=np.int64)
    for r, b in enumerate(fours):
        pa = sorted(set(b) & set(A)); pb = sorted(set(b) & set(B))
        M[r, idx[pa[0]]] += 1; M[r, idx[pa[1]]] += 1; M[r, idx[pb[0]]] -= 1; M[r, idx[pb[1]]] -= 1
    return M, A, B


def kernel_mod_p(M, p):
    M = M.copy() % p; m, n = M.shape
    piv = []; row = 0
    for col in range(n):
        pr = next((r for r in range(row, m) if M[r, col] % p), None)
        if pr is None:
            continue
        M[[row, pr]] = M[[pr, row]]
        inv = pow(int(M[row, col]), -1, p)
        M[row] = (M[row] * inv) % p
        for r in range(m):
            if r != row and M[r, col] % p:
                M[r] = (M[r] - M[r, col] * M[row]) % p
        piv.append(col); row += 1
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(n, dtype=np.int64); v[f] = 1
        for i, c in enumerate(piv):
            v[c] = (-M[i, f]) % p
        basis.append(v)
    sols = []
    for coeffs in itertools.product(range(p), repeat=len(basis)):
        v = np.zeros(n, dtype=np.int64)
        for c, b in zip(coeffs, basis):
            v = (v + c * b) % p
        sols.append(tuple(int(x) for x in v))
    return sols


def torsion_solutions(rec):
    M, A, B = system_matrix(rec)
    k2 = kernel_mod_p(M, 2); k5 = kernel_mod_p(M, 5)
    sols10 = set()
    for s2 in k2:
        for s5 in k5:
            sols10.add(tuple((5 * a + 6 * b) % 10 for a, b in zip(s2, s5)))
    assert all(np.all((M @ np.array(s)) % 10 == 0) for s in sols10)
    good = sorted(s for s in sols10 if len(set(s[:5])) == 5 and len(set(s[5:])) == 5 and s[0] == 0)
    return M, A, B, len(k2), len(k5), len(sols10), good


if __name__ == "__main__":
    data = json.load(open('n10_table_C.json')); rec = data['results'][2]
    M, A, B, n2, n5, n10, good = torsion_solutions(rec)
    print(f"kernel sizes: mod 2: {n2}, mod 5: {n5}; solutions mod 10: {n10} (expected 100 = 10 rotations x 10 torsion classes)")
    print(f"solutions with a_0 = 0, distinct a's and distinct b's: {len(good)}")
    for s in good:
        a, b = s[:5], s[5:]
        print("  angles in units of 2pi/10: a =", a, " b =", b,
              "| a regular pentagon:", sorted(a) in ([0, 2, 4, 6, 8], [1, 3, 5, 7, 9]),
              "| b regular pentagon:", sorted(b) in ([0, 2, 4, 6, 8], [1, 3, 5, 7, 9]),
              "| B rotated by 36 deg relative to A:", (b[0] - a[0]) % 2 == 1)
