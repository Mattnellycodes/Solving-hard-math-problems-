"""Apply the theorem-only necessary conditions (and, separately, the cited-table orchard caps) to
every (F, L) class produced by stage 2.  usage: python3 postfilter.py runs/stage2_*.json"""
import json
import sys
from lib import mask, bits, run_filters, circles_count, check_structure, T3CAP_THEOREM, T3CAP_CITED

survivors = []
tot = 0
for fname in sys.argv[1:]:
    data = json.load(open(fname))
    print(f"== {fname}: {len(data)} (F,L) classes")
    for i, r in enumerate(data):
        F = [mask(b) for b in r["F"]]
        L = [mask(l) for l in r["L"]]
        check_structure(F, L)
        c = circles_count(F, L)
        v = run_filters(F, L, caps=T3CAP_THEOREM)
        vc = run_filters(F, L, caps=T3CAP_CITED)
        kinds = sorted({x[0] for x in v})
        kindsc = sorted({x[0] for x in vc})
        print(f"  #{i}: D={r['D']} l={r['l']} circles={c} b4={sum(1 for b in F if len(bits(b)) == 4)} "
              f"theorem-kills={kinds if kinds else 'NONE (survivor)'}  cited-kills={kindsc if kindsc else 'NONE'}")
        tot += 1
        if not v:
            survivors.append({"file": fname, "index": i, "F": r["F"], "L": r["L"], "circles": c})
print(f"total (F,L) classes: {tot}; survivors of theorem-only filters: {len(survivors)}")
for s in survivors:
    print("SURVIVOR:", s)
json.dump(survivors, open("runs/survivors.json", "w"), indent=1)
