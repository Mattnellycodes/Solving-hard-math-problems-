"""Aut(C-2) and the orbits of its line sets with >= 10 lines (sg caps)."""
import json, mo
C2 = json.load(open('c2_rigidity.json'))
blocks = [mo.mask(b) for b in C2['blocks']]
aut = mo.automorphisms(10, blocks)
print("|Aut(C-2)| =", len(aut))
caps = mo.Caps(10, 'sg')
Ls = mo.line_search(10, blocks, caps, min_size=10)
print("line sets with |L| >= 10:", len(Ls), "sizes:", sorted(len(L) for L in Ls))
canon = {}
for L in Ls:
    key = min(tuple(sorted(mo.apply_perm(s, l) for l in L)) for s in aut)
    canon.setdefault(key, []).append(L)
print("orbits under Aut(C-2):", len(canon))
orbits = []
for key, members in canon.items():
    L = [mo.bits(l) for l in key]
    comp = sorted((len(l) for l in L), reverse=True)
    print(f"  orbit size {len(members)}, |L| = {len(L)}, line sizes {comp}: {L}")
    orbits.append({'size': len(members), 'lines': L})
json.dump({'aut_order': len(aut), 'orbits': orbits}, open('c2_linesets.json', 'w'), indent=1)
