"""Bijection check: the (F,L) classes of my stage 2 (all skeletons) vs the 37 structures stored by
n10-continue in runs/post_all.json, by canonical form of the Möbius structure with inf marked.
Also applies MY filters to THEIR structures and compares kill sets with their recorded kills."""
import json, glob
from lib import mask, bits, canon, mobius_blocks, circles_count, run_filters, check_structure, T3CAP_CITED, T3CAP_THEOREM
mine = {}
for f in sorted(glob.glob("runs/stage2_*.json")):
    if "lb74" in f or "nocap" in f:
        continue
    for i, r in enumerate(json.load(open(f))):
        F = [mask(b) for b in r["F"]]; L = [mask(l) for l in r["L"]]
        cf = canon(mobius_blocks(F, L) + [1 << 10], 11)[0]
        mine.setdefault(cf, []).append((f, i, circles_count(F, L)))
print("my (F,L) classes:", len(mine), " (duplicates across files:", sum(len(v) - 1 for v in mine.values()), ")")
theirs = json.load(open("../n10-continue/runs/post_all.json"))
print("their structures:", len(theirs))
hit = 0; mism = []
for r in theirs:
    F = [mask(b) for b in r["blocks"]]; L = [mask(l) for l in r["lines"]]
    check_structure(F, L)
    cf = canon(mobius_blocks(F, L) + [1 << 10], 11)[0]
    ok = cf in mine
    hit += ok
    v = run_filters(F, L, caps=T3CAP_THEOREM)
    kinds = sorted({x[0] for x in v})
    print(f"  their #{r['index']}: count={r['count']} in mine: {ok}  my theorem-kills={kinds}  their kills={r['theorem_kills']}")
    if not ok:
        mism.append(r["index"])
print("their structures found in my enumeration:", hit, "of", len(theirs), " missing:", mism)
