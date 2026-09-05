"""Integer-coordinate antipodal configurations attaining f(n) = C(n-1,2)+1-floor((n-1)/2), n = 9..16,
certified exactly (two independent counters).  Points: n-1 lattice points on the circle
x^2 + y^2 = 625 in antipodal pairs (plus one unpaired point when n-1 is odd) and the centre."""
import json, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from exact_verify import certify
from circles_exact import count_circles, formula

PAIRS = [((25, 0), (-25, 0)), ((0, 25), (0, -25)), ((7, 24), (-7, -24)), ((24, 7), (-24, -7)),
         ((15, 20), (-15, -20)), ((20, 15), (-20, -15)), ((-7, 24), (7, -24)), ((-24, 7), (24, -7)),
         ((-15, 20), (15, -20)), ((-20, 15), (20, -15))]


def antipodal(n):
    pts = [(0, 0)]
    m = n - 1
    for a, b in PAIRS:
        if len(pts) - 1 + 2 <= m:
            pts += [a, b]
    if len(pts) - 1 < m:
        # one unpaired point: take the first point of the next pair
        pts.append(PAIRS[(m) // 2][0]) if (m % 2 == 1) else None
    assert len(pts) == n
    return pts


if __name__ == "__main__":
    out = {}
    for n in range(9, 17):
        pts = antipodal(n)
        c1 = count_circles(pts)
        r = certify(pts)
        assert c1 == r['circles'] == formula(n), (n, c1, r['circles'], formula(n))
        print(f"n={n}: {c1} circles = f(n) = {formula(n)}; blocks: {r['circles']} circles (sizes {sorted(set(r['circle_sizes']))}), {r['lines']} lines; points {pts}")
        out[n] = dict(points=pts, circles=c1, formula=formula(n), lines=r['lines'])
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs', 'antipodal_certified.json'), 'w'), indent=1)
