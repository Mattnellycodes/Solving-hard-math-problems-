#!/usr/bin/env python3
"""Print, for each candidate structure, the line sets surviving tier A+B (theorems only) and
tier A+B+C (with cited tables), with their circle counts."""
import sys, json
from math import comb
from hcheck import check_structure, cumulative
fn = sys.argv[1]
d = json.load(open(fn)); n = d["n"]
for i, r in enumerate(d["results"]):
    rep = check_structure(n, r["blocks"], r["line_sets"])
    cum = cumulative(rep)
    if not cum['A']:
        continue
    print(f"structure {i}: sizes={r['sizes']} degrees={r['degrees']} D={r['D']}")
    print(f"   blocks={r['blocks']}")
    for t, name in (('B', 'A+B (theorems only)'), ('C', 'A+B+C (with cited tables)')):
        print(f"   surviving line sets, tier {name}: {len(cum[t])}")
        for j in cum[t]:
            L = r["line_sets"][j]
            print(f"      count={comb(n,3)-r['D']-len(L)} lines={L}")
