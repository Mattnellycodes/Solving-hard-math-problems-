"""Post-processing of the candidate structures produced by run_enum.py.
usage: python3 postcheck.py candidates.json [--miquel]
For every candidate family F (count <= target possible):
  * all line sets L with |L| >= need - D(F) are enumerated (mode caps of the run);
  * the derived structure at every point p (lines B\\{p}, B in F, p in B) and the line structure (P, L)
    are tested under two assumption sets:
      A1 = hereditary Sylvester-Gallai (theorem) + no (8_3) (proved in this session; Fano is a special
           case of hereditary SG);
      A2 = A1 + orchard numbers t3(m) on every subset (Burr-Gruenbaum-Sloane 1974, cited);
  * optionally Miquel closure (classical theorem of the real inversive plane).
A structure survives under an assumption set iff its derived structures pass and at least one line set
with |L| >= need - D passes.  The hereditary-SG verdicts are cross-checked against the independent
implementation ../../theory/hereditary_sg.py (sg_violations).
"""
import sys, json, math, itertools, time, os
import mo
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'theory'))
try:
    from hereditary_sg import sg_violations as theory_sg_violations
except Exception as e:
    theory_sg_violations = None
    print("warning: theory/hereditary_sg.py not importable:", e)

path = sys.argv[1]
do_miquel = '--miquel' in sys.argv
data = json.load(open(path))
n, mode, target, need = data['n'], data['mode'], data['target'], data['need']
caps = mo.Caps(n, mode)
cands = data['candidates']
print(f"{path}: n={n} mode={mode} target={target} need={need} candidates={len(cands)} {caps}")
ASSUMPTIONS = {'A1(SG+8_3)': dict(t3=None, check83=True), 'A2(SG+8_3+t3)': dict(t3=mo.T3, check83=True)}
summary = []
for ci, r in enumerate(cands):
    t0 = time.time()
    F = [mo.mask(b) for b in r['blocks']]
    D = mo.deficit(F)
    assert D == r['D']
    min_lines = need - D
    Ls = mo.line_search(n, F, caps, min_size=min_lines)
    ell_max = max((len(L) for L in Ls), default=0)
    print(f"\n=== candidate {ci}: sizes={r['sizes']} D={D} degrees={r['degrees']} ellmax={ell_max} "
          f"count_min={math.comb(n,3)-D-ell_max}; line sets with |L|>={min_lines}: {len(Ls)}")
    print(f"    blocks={r['blocks']}")
    verdict = {}
    for aname, akw in ASSUMPTIONS.items():
        # derived structures
        rep = mo.check_structure(n, F, **akw)
        dead_points = sorted(rep)
        for p in dead_points:
            kind, S, info = rep[p][0]
            print(f"    [{aname}] derived structure at point {p}: {kind} violation on subset {list(S)} {info if kind != 'SG' else '(all pairs covered) lines ' + str(info)}")
        # line structures
        good_L = []
        kinds = {}
        for L in Ls:
            v = mo.hereditary_violations(range(n), [frozenset(mo.bits(l)) for l in L], **akw)
            if v:
                kinds[v[0][0]] = kinds.get(v[0][0], 0) + 1
            else:
                good_L.append(L)
        print(f"    [{aname}] line sets: {len(Ls)} total, {len(good_L)} pass, violations by kind {kinds}")
        if dead_points:
            verdict[aname] = f"EXCLUDED (derived structure at points {dead_points})"
        elif not good_L:
            verdict[aname] = "EXCLUDED (every line set with |L| >= %d violates)" % min_lines
        else:
            best = max(len(L) for L in good_L)
            verdict[aname] = f"SURVIVES (surviving line sets {len(good_L)}, max |L|={best}, min count={math.comb(n,3)-D-best})"
        print(f"    [{aname}] => {verdict[aname]}")
    # cross-check of the hereditary-SG verdicts with the theory agent's implementation
    if theory_sg_violations is not None:
        mine = {p for p in range(n) if any(v[0] == 'SG' for v in mo.hereditary_violations([q for q in range(n) if q != p], mo.derived_lines(n, F, p), t3=None, check83=False, first_only=False))}
        theirs = {p for p in range(n) if theory_sg_violations(set(range(n)) - {p}, mo.derived_lines(n, F, p))}
        print(f"    cross-check hereditary SG (derived, points violating): mine={sorted(mine)} theory/hereditary_sg.py={sorted(theirs)} {'AGREE' if mine == theirs else 'DISAGREE'}")
    miq = None
    if do_miquel:
        mv = mo.miquel_violations(n, F)
        miq = bool(mv)
        if mv:
            v = mv[0]
            print(f"    [Miquel] violation: cube labelling {v['labelling']} faces {v['faces']} -> face {v['missing_face']} {v['faces'][v['missing_face']]} lies in no block  => EXCLUDED by Miquel's theorem")
        else:
            print(f"    [Miquel] Miquel-closed (no violation)")
    summary.append({'index': ci, 'sizes': r['sizes'], 'D': D, 'degrees': r['degrees'], 'ellmax': ell_max,
                    'count_min': math.comb(n, 3) - D - ell_max, 'blocks': r['blocks'], 'verdict': verdict,
                    'miquel_violation': miq, 'n_line_sets': len(Ls), 'time': round(time.time() - t0, 1)})
print("\n==== SUMMARY ====")
for s in summary:
    print(f"candidate {s['index']}: sizes={s['sizes']} D={s['D']} ellmax={s['ellmax']} count_min={s['count_min']} degrees={s['degrees']}")
    for a, v in s['verdict'].items():
        print(f"    {a}: {v}")
    if s['miquel_violation'] is not None:
        print(f"    Miquel: {'VIOLATION (excluded)' if s['miquel_violation'] else 'closed'}")
for aname in ASSUMPTIONS:
    surv = [s['index'] for s in summary if s['verdict'][aname].startswith('SURVIVES')]
    print(f"survivors under {aname}: {surv}")
if do_miquel:
    for aname in ASSUMPTIONS:
        surv = [s['index'] for s in summary if s['verdict'][aname].startswith('SURVIVES') and not s['miquel_violation']]
        print(f"survivors under {aname} + Miquel: {surv}")
json.dump(summary, open(path.replace('.json', '_postcheck.json'), 'w'), indent=1)
