#!/usr/bin/env python3
# Hand-derived rigorous bounds for pi and gamma with exact rationals (no library constants):
#  pi via Machin's formula (alternating series, error <= first omitted term);
#  gamma via Euler-Maclaurin with the signed remainder of Lemma 4 (theta in (0,1)), N = 1024, m = 4, ln 2 = 2 atanh(1/3).
from fractions import Fraction as Fr
def atan_bounds(x, N):   # x = Fr, 0<x<=1: partial sums of alternating series with decreasing terms
    s = Fr(0)
    for j in range(N): s += (-1)**j * x**(2*j+1) / (2*j+1)
    t = x**(2*N+1) / (2*N+1)
    return s - t, s + t
a5 = atan_bounds(Fr(1,5), 40); a239 = atan_bounds(Fr(1,239), 20)
pi_lo = 16*a5[0] - 4*a239[1]; pi_hi = 16*a5[1] - 4*a239[0]
print("pi in [%.30f, %.30f]" % (pi_lo, pi_hi)); print("pi <= 3.1415926536:", pi_hi <= Fr('3.1415926536'), "; pi >= 3.1415926535:", pi_lo >= Fr('3.1415926535'))
K = 45
S = sum(Fr(1, (2*k+1)*3**(2*k+1)) for k in range(K))
ln2_lo = 2*S; ln2_hi = 2*S + 2*Fr(9,8)*Fr(1, (2*K+1)*3**(2*K+1))
N = 1024; H = sum(Fr(1, i) for i in range(1, N+1))
B = {2: Fr(1,6), 4: Fr(-1,30), 6: Fr(1,42), 8: Fr(-1,30), 10: Fr(5,66)}
X = lambda ln: H - ln - Fr(1, 2*N) + sum(B[2*j]/(2*j*N**(2*j)) for j in range(1, 5))
# gamma = X(ln N) + theta*B_10/(10 N^10), theta in (0,1)  [from (3.2) with d=N, m=4]
g_lo = X(10*ln2_hi); g_hi = X(10*ln2_lo) + B[10]/(10*N**10)
print("gamma in [%.30f, %.30f]" % (g_lo, g_hi)); print("gamma >= 0.5772156649:", g_lo >= Fr('0.5772156649'), "; gamma <= 0.5772156650:", g_hi <= Fr('0.5772156650'))
print("width pi:", float(pi_hi-pi_lo), " width gamma:", float(g_hi-g_lo))
