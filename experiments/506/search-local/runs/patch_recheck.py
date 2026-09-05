p = '/home/user/Solving-hard-math-problems-/experiments/506/search-local/recheck_results.py'
src = open(p).read()
old = """    if tc < r['proxy']:
        ok = False
        try:
            ex = RS.exact_coords(U, S)
            if ex is not None:"""
new = """    if tc < r['proxy']:
        ok = False
        try:
            ex = RS.exact_coords(U, S)
            if tc > formula(r['n']) + 2:
                ex = None; note += "; far above f(n): certification skipped"
            if ex is not None:"""
assert old in src
open(p, 'w').write(src.replace(old, new))
print("patched")
