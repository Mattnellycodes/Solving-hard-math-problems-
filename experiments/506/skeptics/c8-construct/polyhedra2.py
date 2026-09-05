import numpy as np, itertools, math
import polyhedra as ph
from robust import evaluate, normalise
rng = np.random.default_rng(7)
def sphere_proj(V, view):
    V = np.array(V, float); V /= np.linalg.norm(V, axis=1)[:, None]
    N = np.array(view, float); N /= np.linalg.norm(N)
    if max(V @ N) > 0.9: return None   # view too close to a vertex -> huge coordinates
    z = N; xx = np.cross(z, [0.3, 0.5, 0.8]); xx /= np.linalg.norm(xx); yy = np.cross(z, xx)
    W = V @ np.array([xx, yy, z]).T
    return np.array([[w[0]/(1-w[2]), w[1]/(1-w[2])] for w in W])
polys = {}
for h in [0.3, 0.5, 0.7071, 1.0, 1.5, 1.7320508]:
    polys[f"square antiprism h={h}"] = [(math.cos(k*math.pi/2), math.sin(k*math.pi/2), h) for k in range(4)] + [(math.cos(k*math.pi/2+math.pi/4), math.sin(k*math.pi/2+math.pi/4), -h) for k in range(4)]
    polys[f"square prism h={h}"] = [(math.cos(k*math.pi/2), math.sin(k*math.pi/2), h) for k in range(4)] + [(math.cos(k*math.pi/2), math.sin(k*math.pi/2), -h) for k in range(4)]
    polys[f"triangular prism+2 apex h={h}"] = [(math.cos(2*k*math.pi/3), math.sin(2*k*math.pi/3), h) for k in range(3)] + [(math.cos(2*k*math.pi/3), math.sin(2*k*math.pi/3), -h) for k in range(3)] + [(0,0,1),(0,0,-1)]
    polys[f"hexagon+2 poles h={h}"] = [(math.cos(k*math.pi/3), math.sin(k*math.pi/3), h) for k in range(6)] + [(0,0,1),(0,0,-1)]
    polys[f"twisted trig prism+2 h={h}"] = [(math.cos(2*k*math.pi/3), math.sin(2*k*math.pi/3), h) for k in range(3)] + [(math.cos(2*k*math.pi/3+math.pi/3), math.sin(2*k*math.pi/3+math.pi/3), -h) for k in range(3)] + [(0,0,1),(0,0,-1)]
    polys[f"rect prism 1x2 h={h}"] = [(x, y, h) for x in (-1, 1) for y in (-2, 2)] + [(x, y, -h) for x in (-1, 1) for y in (-2, 2)]
    polys[f"rect prism 1x1.7320508 h={h}"] = [(x, y, h) for x in (-1, 1) for y in (-1.7320508, 1.7320508)] + [(x, y, -h) for x in (-1, 1) for y in (-1.7320508, 1.7320508)]
polys["cube"] = list(itertools.product([-1, 1], repeat=3))
polys["cube 2 opposite vertices replaced by face centres?"] = [(-1,-1,-1),(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1),(1,1,-1),(1,-1,1),(-1,1,1)]
summary = {}
for name, V in polys.items():
    b = 99; nb = None; unstable = 0
    for trial in range(60):
        P = sphere_proj(V, rng.normal(size=3))
        if P is None: continue
        r = evaluate(P)
        if r is None: continue
        if r[0] == "UNSTABLE": unstable += 1; continue
        b = min(b, r[0]); nb = r[2]
    print(f"{name}: |B|={nb} best over views+inversions = {b}  (unstable views: {unstable})")
