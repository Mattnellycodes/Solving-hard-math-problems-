"""Rational (integer) antipodal configurations attaining f(n) for n = 9..16, certified exactly:
n-1 points on the unit circle (k antipodal pairs, plus one extra rational point when n-1 is odd)
and the centre.  Output: integer coordinates."""
import sys, math, json
from fractions import Fraction as Fr
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506')
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/search-local')
from circles_exact import count_circles, formula, all_on_one_circle_or_line
from verify_exact import count_rational

def rat_circle_points(k):
    """k rational points on the unit circle with pairwise distinct antipodal pairs"""
    pts = []
    t = 1
    while len(pts) < k:
        x = Fr(1 - t * t, 1 + t * t); y = Fr(2 * t, 1 + t * t)
        if (x, y) not in pts and (-x, -y) not in pts:
            pts.append((x, y))
        t += 1
    return pts

out = {}
for n in range(9, 17):
    m = n - 1
    k = m // 2
    base = rat_circle_points(k + 1)
    pts = []
    for (x, y) in base[:k]:
        pts += [(x, y), (-x, -y)]
    if m % 2 == 1:
        pts.append(base[k])
    pts.append((Fr(0), Fr(0)))
    assert len(pts) == n and len(set(pts)) == n
    c1 = count_circles(pts)
    c2, sizes, lsizes = count_rational(pts)
    assert c1 == c2 == formula(n), (n, c1, c2)
    assert not all_on_one_circle_or_line(pts)
    den = 1
    for x, y in pts:
        den = den * x.denominator // math.gcd(den, x.denominator)
        den = den * y.denominator // math.gcd(den, y.denominator)
    coords = [(int(x * den), int(y * den)) for x, y in pts]
    out[n] = dict(n=n, circles=c1, formula=formula(n), coords=coords, circle_sizes=sizes, line_sizes=lsizes)
    print(f"n={n:2d} f={formula(n):3d} circles={c1:3d} sizes={sizes} lines={lsizes} coords={coords}")
json.dump(out, open('/home/user/Solving-hard-math-problems-/experiments/506/search-local/runs/antipodal_rational.json', 'w'), indent=1)
