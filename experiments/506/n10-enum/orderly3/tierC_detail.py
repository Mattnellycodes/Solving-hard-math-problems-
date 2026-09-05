#!/usr/bin/env python3
"""For the line sets surviving tiers A+B+K, list all tier-C violations (which cited value kills them)."""
import json
from math import comb
from hcheck import check_structure, cumulative
from oe import list_to_mask, mask_to_list
d = json.load(open('n10_sg_t32.json')); n = d['n']
for i in (4, 5, 6, 7, 8):
    r = d['results'][i]
    rep = check_structure(n, r['blocks'], r['line_sets'])
    cum = cumulative(rep)
    print(f"structure {i}: sizes={r['sizes']} degrees={r['degrees']}: A+B+K survivors {cum['K']}, A+B+K+C survivors {cum['C']}")
    for j in cum['K']:
        L = r['line_sets'][j]
        print(f"   line set {j} (count {comb(n,3)-r['D']-len(L)}): lines={L}")
        lv = rep['line_sets'][j]['C']
        if not lv:
            print("      passes tier C")
        for (q, ex, k) in lv:
            print(f"      tier C violation at point {q} ({k} violations), e.g. {ex}")
