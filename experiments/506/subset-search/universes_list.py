"""The universes used by the driver.  Each entry: (name, constructor thunk)."""
import sympy as sp
import universe as uv
s2, s3 = sp.sqrt(2), sp.sqrt(3)
phi = (1 + sp.sqrt(5)) / 2

def sec(m): return 1 / sp.cos(sp.pi / m)

UNIVERSES = {
    # integer grids
    'grid4': lambda: uv.grid(4),
    'grid5': lambda: uv.grid(5),
    'grid6': lambda: uv.grid(6),
    'grid7': lambda: uv.grid(7),
    'grid9': lambda: uv.grid(9),
    'grid5half': lambda: uv.grid_half(5),
    'grid7half': lambda: uv.grid_half(7),
    'grid4x8': lambda: uv.grid(4, 8),
    'grid3x12': lambda: uv.grid(3, 12),
    # lattice points on concentric circles about the origin (many rational points per circle)
    'latcirc_a': lambda: uv.lattice_circles([25, 50, 100, 200, 400]),
    'latcirc_b': lambda: uv.lattice_circles([65, 130, 260, 325, 650]),
    'latcirc_c': lambda: uv.lattice_circles([25, 65, 100, 125, 169, 325, 425]),
    'latcirc_d': lambda: uv.lattice_circles([1, 2, 4, 5, 8, 9, 10, 13, 16, 18, 20, 25, 26, 29, 32, 34, 36, 37, 40, 41, 45, 49, 50]),
    # concentric regular polygons (both rotations) + centre
    'poly3': lambda: uv.polygons(3, [1, 2, 3, 4, sp.Rational(1, 2), s3, 2 * s3, 3 * s3]),
    'poly4': lambda: uv.polygons(4, [1, 2, 3, 4, s2, 2 * s2, 3 * s2, 1 + s2]),
    'poly5': lambda: uv.polygons(5, [1, 2, 3, phi, phi ** 2, sec(5), sec(5) ** 2, sec(10)]),
    'poly6': lambda: uv.polygons(6, [1, 2, 3, 4, s3, 2 * s3, sec(6), 2 * sec(6)]),
    'poly8': lambda: uv.polygons(8, [1, 2, 3, s2, 1 + s2, sec(8), sec(8) ** 2, 2 * sec(8)]),
    'poly10': lambda: uv.polygons(10, [1, 2, phi, sec(10), phi * sec(10), sec(5)]),
    'poly12': lambda: uv.polygons(12, [1, 2, s3, sec(12), 2 * sec(12), sec(6)]),
    'poly7': lambda: uv.polygons(7, [1, 2, sec(7), sec(7) ** 2, sec(14), 2 * sec(7)]),
    'poly9': lambda: uv.polygons(9, [1, 2, sec(9), sec(9) ** 2, sec(18), 2 * sec(9)]),
    # points on a circle plus lattice points
    'c24grid2': lambda: uv.circle_plus_grid(24, 1, 2),
    'c12grid3': lambda: uv.circle_plus_grid(12, 1, 3),
    'c16grid2': lambda: uv.circle_plus_grid(16, 1, 2),
    'c24R2grid3': lambda: uv.circle_plus_grid(24, 2, 3),
    'c20grid2': lambda: uv.circle_plus_grid(20, 1, 2),
    'c8grid2': lambda: uv.circle_plus_grid(8, 1, 2),
    'c10grid2': lambda: uv.circle_plus_grid(10, 1, 2),
}
