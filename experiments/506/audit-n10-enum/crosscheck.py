"""Brute-force S_n canonical forms: compare my candidate classes with the other agent's JSON."""
import sys, json
from mob import PermTable, canon_full
mine = json.load(open(sys.argv[1])); theirs = json.load(open(sys.argv[2]))
n = mine['n']; PT = PermTable(n)
cm = {}
for i, c in enumerate(mine['results']): cm.setdefault(canon_full(PT, c['blocks']), []).append(('mine', i, c['count_min'], len(c['line_sets'])))
ct = {}
for i, c in enumerate(theirs['results']): ct.setdefault(canon_full(PT, c['blocks']), []).append(('theirs', i, c['count_min'], len(c['line_sets'])))
print(f"my classes: {len(cm)} (from {len(mine['results'])} structures); their classes: {len(ct)} (from {len(theirs['results'])})")
print("in both:", sum(1 for k in cm if k in ct), " only mine:", [cm[k] for k in cm if k not in ct], " only theirs:", [ct[k] for k in ct if k not in cm])
for k in cm:
    if k in ct: print("  match:", cm[k], "<->", ct[k])
