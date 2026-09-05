"""Recount the 8 antipodal records (coordinates copied from the search-inversion report text)."""
from count506 import analyse, f
base = [(0,0),(25,0),(-25,0),(0,25),(0,-25),(7,24),(-7,-24),(24,7),(-24,-7),(15,20),(-15,-20),(20,15),(-20,-15),(-7,24),(7,-24),(-24,7)]
ok = True
for n in range(9, 17):
    pts = base[:n]
    r = analyse(pts)
    circ_pts = [p for p in pts if p[0]**2 + p[1]**2 == 625]
    pairs = sum(1 for p in circ_pts if (-p[0], -p[1]) in circ_pts) // 2
    good = (r['circles'] == f(n)) and not r['degenerate'] and r['check_formula'] == r['circles']
    ok &= good
    print(f"n={n:2d}: circles={r['circles']:3d} f(n)={f(n):3d} lines={r['lines']} line_sizes={r['line_sizes']} "
          f"circle_sizes: max={r['circle_sizes'][0]}, #3s={r['circle_sizes'].count(3)} | on circle: {len(circ_pts)} "
          f"antipodal pairs: {pairs} | degenerate={r['degenerate']} | OK={good}")
print("ALL 8 RECORDS EQUAL f(n) AND NON-DEGENERATE:", ok)
# generic-structure check: exactly one (n-1)-circle, floor((n-1)/2) 3-point lines, all other circles 3-point
for n in range(9, 17):
    r = analyse(base[:n])
    assert r['circle_sizes'][0] == n - 1 and r['circle_sizes'][1:] == [3] * (r['circles'] - 1)
    assert r['line_sizes'] == [3] * ((n - 1) // 2)
print("structure: one (n-1)-circle + floor((n-1)/2) three-point lines + 3-point circles only: verified for n=9..16")
