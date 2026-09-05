"""Corrected verdicts for c2_exact.json (sympy intervals are ((a, b), multiplicity)).
A (type, line set) case is EXCLUDED iff (a) the gcd of the g_k, after removing the trivial factors rho
and rho^2 - 1, is constant (so for rho not in {0, 1} the blocks of L never all pass through the radical
centre), and (b) Delta has no positive real root other than rho = 1 (so the radical centre is always
defined).  Root intervals are exact isolating intervals of the rational norm polynomials."""
import json, sys
from fractions import Fraction as Fr
rep = []
for fn in (sys.argv[1:] or ['c2_exact.json']):
    rep += json.load(open(fn))
def parse_iv(t):
    a, b = t
    return (Fr(a.strip('()') .split(',')[0]) if '(' in a else Fr(a), Fr(b))
bad = []
for e in rep:
    droots = e['delta_norm_positive_real_roots'] or []
    # each entry is (str(interval), str(multiplicity)); interval printed as "(a, b)"
    others = []
    for ivs, mult in droots:
        a, b = [Fr(x) for x in ivs.strip('()').split(',')]
        if not (a == 1 and b == 1):
            others.append((a, b))
    gred = e['reduced_gcd_degree']
    ok = (gred == 0) and not others
    if gred > 0:
        rr = e.get('reduced_gcd_norm_positive_real_roots', [])
        rr2 = []
        for ivs, mult in rr:
            a, b = [Fr(x) for x in ivs.strip('()').split(',')]
            if not (a == 1 and b == 1):
                rr2.append((a, b))
        ok = (not rr2) and not others
    e['verdict_corrected'] = 'EXCLUDED' if ok else 'NEEDS ANALYSIS'
    if not ok:
        bad.append(e)
types = sorted({e['type_offset'] for e in rep}); lsets = sorted({e['lineset'] for e in rep})
print(f"cases: {len(rep)} (types {types}, line sets {len(lsets)}); EXCLUDED: {len(rep) - len(bad)}; needing analysis: {len(bad)}")
for e in bad:
    print("  NEEDS ANALYSIS:", {k: e[k] for k in ('type_offset', 'lineset', 'reduced_gcd_degree', 'reduced_gcd', 'delta_norm_positive_real_roots')})
print("complete (2 types x 31 line sets)?", len(rep) == 62)
if len(rep) == 62 and not bad:
    print("RESULT: no realisation of C-2 admits a point outside P lying on >= 10 blocks; every planar point set with block structure C-2 determines >= 42 - 9 = 33 circles.")
