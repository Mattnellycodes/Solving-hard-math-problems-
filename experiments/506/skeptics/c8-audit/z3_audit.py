"""Independent z3 check: is there an abstract structure (F, L) on 8 points with <= TARGET circles
under (C1)-(C3)?  Same relaxation as cpsat_audit.py but a different solver and encoding (PB).
Usage: python3 z3_audit.py MODE TARGET"""
import sys, itertools, time
from math import comb
import z3
N = 8
OVALS = {"strong": {7: 3, 8: 4}, "weak": {7: 1, 8: 1}, "none": {7: 0, 8: 0}}
mode, target = sys.argv[1], int(sys.argv[2])
o = OVALS[mode]; capd, capl = comb(N - 1, 2) - o[N - 1], comb(N, 2) - o[N]
subs = lambda k: [frozenset(c) for c in itertools.combinations(range(N), k)]
big = [B for k in range(4, N) for B in subs(k)]; tri = subs(3)
x = {B: z3.Bool("x%s" % "".join(map(str, sorted(B)))) for B in big}
z = {t: z3.Bool("z%s" % "".join(map(str, sorted(t)))) for t in tri}
y = {S: z3.Bool("y%s" % "".join(map(str, sorted(S)))) for S in big + tri}
s = z3.Solver()
for t in tri:   # exactly one block contains t
    s.add(z3.PbEq([(z[t], 1)] + [(x[B], 1) for B in big if t <= B], 1))
for B1, B2 in itertools.combinations(big, 2):
    if len(B1 & B2) >= 3: s.add(z3.Or(z3.Not(x[B1]), z3.Not(x[B2])))
for S in big: s.add(z3.Implies(y[S], x[S]))
for t in tri: s.add(z3.Implies(y[t], z[t]))
for S1, S2 in itertools.combinations(big + tri, 2):
    if len(S1 & S2) >= 2: s.add(z3.Or(z3.Not(y[S1]), z3.Not(y[S2])))
for p in range(N):
    s.add(z3.PbLe([(x[B], comb(len(B) - 1, 2)) for B in big if p in B], capd))
s.add(z3.PbLe([(y[S], comb(len(S), 2)) for S in big + tri], capl))
# circles = #x + #z - #y <= target  <=>  #x + #z + #(not y) <= target + #S
allS = big + tri
s.add(z3.PbLe([(x[B], 1) for B in big] + [(z[t], 1) for t in tri] + [(z3.Not(y[S]), 1) for S in allS],
              target + len(allS)))
t0 = time.time(); r = s.check()
print(f"mode={mode} target={target} (capd={capd}, capl={capl}): {r}  [{time.time()-t0:.1f}s]")
if r == z3.sat:
    M = s.model()
    F = sorted(sorted(B) for B in big if z3.is_true(M.eval(x[B], model_completion=True)))
    L = sorted(sorted(S) for S in allS if z3.is_true(M.eval(y[S], model_completion=True)))
    Z = sum(1 for t in tri if z3.is_true(M.eval(z[t], model_completion=True)))
    print("  blocks:", F); print("  lines:", L); print("  #3-blocks:", Z, " circles =", len(F) + Z - len(L))
