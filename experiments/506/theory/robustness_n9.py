"""Re-run the n=9 (target 24) enumeration with weaker ordinary-line bounds and apply the hereditary
SG filter, to show the conclusion c(9) >= 25 does not depend on the exact small values of o(m).
Variant CS: o(m) = ceil(6m/13) for m >= 8 (Csima-Sawyer only), o(7) = 3 (Kelly-Moser).
Variant SG: o(m) = 1 for all m (Sylvester-Gallai only)."""
import json, sys
import mobius_enum as me, mobius_enum2 as me2
from hereditary_sg import check_structure
variants = {"CS": lambda m: {3:3,4:3,5:4,6:3,7:3}.get(m, -(-6*m//13)) if m >= 3 else 0,
            "SG": lambda m: 1 if m >= 3 else 0}
for name, fn in variants.items():
    E = me2.Enumerator2(9, 24, [5,6,7,8], verbose=False, o_func=fn)
    print(f"variant {name}: cap_derived={E.cap_derived} cap_lines={E.cap_lines} ell_max={E.ell_max} D_min={E.D_min}", flush=True)
    res = E.run(include_empty=False)
    print(f"  nodes {E.nodes}, candidates with count<=24: {len(res)}")
    for r in res:
        rep, lrep = check_structure(r)
        killed = sorted(rep)
        # also check for a point of derived degree 8 -> (8_3) (Lemma 2.4): count 4-blocks through p
        print(f"   sizes={r['sizes']} D={r['D']} ell={r['ell_max']} count={r['count_min']} degrees={r['degrees']} "
              f"hereditary-SG violated at points {killed}")
