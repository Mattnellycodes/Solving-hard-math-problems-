#!/usr/bin/env python3
"""Compare candidate lists (JSON files with 'results': [{'blocks': ...}]) up to isomorphism, using the
brute-force S_n canonical form (all n! permutations) -- independent of the orderly machinery."""
import sys, json, time
import numpy as np
from oe import all_perms, canon_bruteforce, list_to_mask

def load(fn):
    d = json.load(open(fn))
    if isinstance(d, dict):
        res = d.get("results", d.get("candidates"))
        res = [{"blocks": r["blocks"], "count_min": r.get("count_min", r.get("count")),
                "ell_max": r.get("ell_max", r.get("ellmax"))} for r in res]
        return d.get("n", 10), res
    return 10, d

def main():
    files = sys.argv[1:]
    n = None
    forms = {}
    perms = None
    for fn in files:
        nn, res = load(fn)
        if perms is None:
            n = nn; perms = all_perms(n)
        forms[fn] = []
        for r in res:
            fam = tuple(sorted(list_to_mask(b) for b in r["blocks"]))
            cf = canon_bruteforce(n, fam, perms)
            forms[fn].append((cf, r.get("count_min"), r.get("ell_max")))
        print(f"{fn}: {len(res)} structures, {len(set(f[0] for f in forms[fn]))} distinct classes", flush=True)
    ref = files[0]
    for fn in files[1:]:
        mine = {f[0]: i for i, f in enumerate(forms[ref])}
        theirs = {f[0]: i for i, f in enumerate(forms[fn])}
        for cf, i in mine.items():
            j = theirs.get(cf)
            print(f"  {ref}[{i}] (count {forms[ref][i][1]}, ell {forms[ref][i][2]}) -> "
                  f"{fn}[{j}]" + (f" (count {forms[fn][j][1]}, ell {forms[fn][j][2]})" if j is not None else "  MISSING"))
        for cf, j in theirs.items():
            if cf not in mine:
                print(f"  {fn}[{j}] has NO counterpart in {ref}")
        same = set(mine) == set(theirs)
        print(f"  => classes {'IDENTICAL' if same else 'DIFFER'} ({len(mine)} vs {len(theirs)})")

if __name__ == '__main__':
    main()
