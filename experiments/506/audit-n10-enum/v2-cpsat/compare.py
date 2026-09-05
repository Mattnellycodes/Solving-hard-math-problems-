"""compare.py mine.json [mine2.json ...] -- theirs.json [theirs2.json ...]
Isomorphism comparison (own canon codes) of structure lists; also prints per-class counts."""
import sys, json
from canon import canon
args = sys.argv[1:]
sep = args.index('--')
mine_files, their_files = args[:sep], args[sep + 1:]
def load(files, key):
    out = {}
    for fn in files:
        d = json.load(open(fn))
        n = d['n']
        for i, r in enumerate(d['results']):
            code = canon(n, [frozenset(B) for B in r['blocks']])[0]
            cnt = r.get('count', r.get('count_min'))
            nls = r.get('n_line_sets', len(r.get('line_sets', [])))
            out.setdefault(code, []).append((fn.split('/')[-1], i, cnt, nls, r.get('survivors')))
    return out
M = load(mine_files, 'mine'); T = load(their_files, 'theirs')
print(f'my classes: {len(M)}; their classes: {len(T)}; common: {len(set(M) & set(T))}')
for c in sorted(set(M) | set(T)):
    print('  BOTH ' if c in M and c in T else ('  MINE ' if c in M else '  THEIRS'), M.get(c), '<->', T.get(c))
