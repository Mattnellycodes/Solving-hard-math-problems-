"""Simulated annealing over a large rational universe: choose k points Q, minimise
min_q (|B(Q)| - deg(q)) (best Euclidean count for n=k-1 after inversion) and also |B(Q)| - #lines (n=k).
Universe: 7x7 grid + rational points on two circles + some half-integer points."""
import sys, random, itertools, time
from fractions import Fraction as Fr
from gridsearch import block_key
def make_universe():
    pts = [(Fr(x), Fr(y)) for x in range(7) for y in range(7)]
    cx, cy = Fr(3), Fr(3)
    for R2, offs in [(Fr(25), [(5,0),(0,5),(3,4),(4,3)]), (Fr(25,4), [(Fr(5,2),0),(0,Fr(5,2)),(Fr(3,2),2),(2,Fr(3,2))])]:
        for dx, dy in offs:
            for sx in (1,-1):
                for sy in (1,-1):
                    p = (cx + sx*dx, cy + sy*dy)
                    if p not in pts: pts.append(p)
    for x in range(7):
        for y in range(7):
            p = (Fr(x) + Fr(1,2), Fr(y) + Fr(1,2))
            if p not in pts: pts.append(p)
    return pts
U = make_universe()
m = len(U)
cache = {}
def key(i, j, k):
    t = (i, j, k) if i < j < k else tuple(sorted((i, j, k)))
    if t not in cache: cache[t] = block_key(U[t[0]], U[t[1]], U[t[2]])
    return cache[t]
def evaluate(Q):
    blocks = {}
    for i, j, k in itertools.combinations(Q, 3):
        blocks.setdefault(key(i, j, k), set()).update((i, j, k))
    nb = len(blocks)
    if nb == 1: return 10**9, 10**9
    nl = sum(1 for kk in blocks if kk[0] == 'L')
    best = 10**9
    for q in Q:
        if any(len(s) == len(Q) - 1 and q not in s for s in blocks.values()): continue
        deg = sum(1 for s in blocks.values() if q in s)
        best = min(best, nb - deg)
    return best, nb - nl
def sa(k, iters, seed):
    rng = random.Random(seed)
    Q = rng.sample(range(m), k)
    cur = evaluate(Q); cur_score = min(cur[0] * 1000 + cur[1], 10**12)
    best = (cur, list(Q))
    T = 2.0
    for it in range(iters):
        i = rng.randrange(k); new = rng.randrange(m)
        if new in Q: continue
        Q2 = Q[:]; Q2[i] = new
        v = evaluate(Q2); sc = min(v[0] * 1000 + v[1], 10**12)
        if sc <= cur_score or rng.random() < pow(2.718, -(sc - cur_score) / (1000 * T)):
            Q, cur, cur_score = Q2, v, sc
            if v[0] < best[0][0] or (v[0] == best[0][0] and v[1] < best[0][1]):
                best = (v, list(Q))
        T = max(0.05, T * 0.9995)
    return best
if __name__ == '__main__':
    k = int(sys.argv[1]); runs = int(sys.argv[2]); iters = int(sys.argv[3])
    print('universe size', m, 'k', k, flush=True)
    overall = None
    for r in range(runs):
        b = sa(k, iters, r)
        print(f'run {r}: best (inversion n={k-1}: {b[0][0]}, direct n={k}: {b[0][1]})  Q={[(str(U[i][0]), str(U[i][1])) for i in b[1]]}', flush=True)
        if overall is None or b[0] < overall[0]: overall = b
    print('OVERALL', overall[0], flush=True)
