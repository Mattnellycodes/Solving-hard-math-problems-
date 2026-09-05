"""run.py n count_target [--skel-min-size k] [--only i,j] : big-block case (skeleton = blocks of size >= 5).
For every skeleton class, CP-SAT enumerates all completions (F4, L) with D + ell >= C(n,3) - count_target,
then structures are classified up to isomorphism with canon(); line sets are re-enumerated independently
(line_sets) and the hereditary checks are run on every (F, L)."""
import sys, time, json, itertools
from math import comb
from pipe import *

n = int(sys.argv[1]); ctarget = int(sys.argv[2])
only = None
if '--only' in sys.argv: only = set(int(t) for t in sys.argv[sys.argv.index('--only') + 1].split(','))
tlimit = 2400
if '--tlimit' in sys.argv: tlimit = int(sys.argv[sys.argv.index('--tlimit') + 1])
caps = caps_for(n, 'sg')
target = comb(n, 3) - ctarget
print(f'n={n} count<={ctarget}: need D+ell>={target}; caps {caps}', flush=True)
t0 = time.time()
skels = enumerate_skeletons(n, caps['capD'])
skels = [s for s in skels if s]     # nonempty (all-4-block case handled by all4.py)
print(f'{len(skels)} nonempty skeleton classes: {sorted(sorted(len(B) for B in s) for s in skels)} [{time.time()-t0:.1f}s]', flush=True)
classes = {}
incomplete = []
for i, skel in enumerate(skels):
    if only is not None and i not in only: continue
    t1 = time.time()
    st, sols = complete(n, skel, target, caps, time_limit=tlimit)
    F4s = {F4 for F4, L in sols}
    print(f'skeleton {i} sizes {sorted(len(B) for B in skel)}: status {st}, {len(sols)} (F4,L) solutions, {len(F4s)} distinct F4 [{time.time()-t1:.1f}s]', flush=True)
    if st not in ('OPTIMAL', 'FEASIBLE') and st != 'INFEASIBLE':
        incomplete.append(i)
    if st == 'FEASIBLE':   # time limit hit while enumerating
        incomplete.append(i)
    for F4 in F4s:
        F = list(skel) + sorted(F4, key=sorted)
        code = canon(n, F)[0]
        if code not in classes:
            classes[code] = dict(skeleton=i, blocks=F, nlab=0, cp_linesets=set())
        classes[code]['nlab'] += 1
        classes[code]['cp_linesets'].update(L for F4b, L in sols if F4b == F4)
print(f'\n{len(classes)} isomorphism classes of rich-block structures with some completion; incomplete skeletons: {incomplete}', flush=True)
results = []
for k, (code, c) in enumerate(sorted(classes.items(), key=lambda kv: kv[0])):
    F = c['blocks']
    D = sum(comb(len(B), 3) - 1 for B in F)
    Ls = line_sets(n, F, caps['capL'], target - D)
    ellmax = max(len(L) for L in Ls) if Ls else -1
    deg = sorted(sum(1 for B in F if p in B) for p in range(n))
    auto = aut_group_order(n, F)[0]
    cnt = comb(n, 3) - D - ellmax
    # cross-check: CP-SAT line sets for the *representative* labelling? (they are for other labellings;
    # compare only counts of line sets per structure via orbit sizes later). Here: max size agrees?
    cpmax = max((len(L) for L in c['cp_linesets']), default=-1)
    rec = dict(idx=k, skeleton=c['skeleton'], sizes=sorted((len(B) for B in F), reverse=True), degrees=deg, aut=auto,
               D=D, ell_max=ellmax, count=cnt, n_line_sets=len(Ls), cp_ellmax=cpmax,
               blocks=[sorted(B) for B in F], line_sets=[sorted(sorted(l) for l in L) for L in Ls])
    # hereditary checks on every line set
    surv = {'A': 0, 'AB': 0, 'ABK': 0, 'ABKC': 0}
    fails_summary = []
    for L in Ls:
        f = full_check(n, F, L)
        okA = not f['A']; okB = okA and not f['B']; okK = okB and not f['K']; okC = okK and not f['C']
        surv['A'] += okA; surv['AB'] += okB; surv['ABK'] += okK; surv['ABKC'] += okC
        fails_summary.append({t: [(str(w[0]), str(w[1])[:120]) for w in f[t]][:1] for t in f})
    rec['survivors'] = surv
    rec['fails'] = fails_summary
    results.append(rec)
    print(f"class {k}: sizes {rec['sizes']} deg {deg} |Aut| {auto} D={D} ell_max={ellmax} (cp {cpmax}) count={cnt} "
          f"#linesets={len(Ls)} survivors(A/AB/ABK/ABKC)={surv} labelled={c['nlab']}", flush=True)
    if cnt <= ctarget:
        print('   blocks', [sorted(B) for B in F])
json.dump(dict(n=n, ctarget=ctarget, incomplete=incomplete, results=results), open(f'out_n{n}_t{ctarget}.json', 'w'), indent=1)
print(f'total time {time.time()-t0:.1f}s')
