"""Compare my candidate list with the earlier agent's list (../n10_sg.json, results[*].blocks) up to
isomorphism, in both directions."""
import sys, json
import mo
mine = json.load(open(sys.argv[1]))['candidates']
theirs = json.load(open(sys.argv[2]))['results']
n = 10
M = [[mo.mask(b) for b in r['blocks']] for r in mine]
T = [[mo.mask(b) for b in r['blocks']] for r in theirs]
print(f"mine: {len(M)} candidates; theirs: {len(T)} candidates")
match_mt = {}
for i, F in enumerate(M):
    js = [j for j, G in enumerate(T) if mo.is_isomorphic(n, F, G)]
    match_mt[i] = js
    print(f"  mine[{i}] sizes={mine[i]['sizes']} D={mine[i]['D']} ellmax={mine[i]['ellmax']} -> theirs {js}"
          + (f" (their ell_max={theirs[js[0]]['ell_max']}, count_min={theirs[js[0]]['count_min']})" if js else " UNMATCHED"))
for j, G in enumerate(T):
    is_ = [i for i, F in enumerate(M) if mo.is_isomorphic(n, F, G)]
    print(f"  theirs[{j}] sizes={theirs[j]['sizes']} D={theirs[j]['D']} ell_max={theirs[j]['ell_max']} -> mine {is_}" + ("" if is_ else " UNMATCHED"))
# pairwise non-isomorphic among mine?
dups = [(i, j) for i in range(len(M)) for j in range(i+1, len(M)) if mo.is_isomorphic(n, M[i], M[j])]
print("isomorphic pairs among my candidates:", dups)
