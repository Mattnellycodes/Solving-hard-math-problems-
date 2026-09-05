"""Incremental post-processing of runs/sa_results.jsonl -> runs/sa_results_checked.jsonl.
Every record whose float 'true' value is below the exact combinatorial proxy is re-evaluated with
the tightened true_count; if still below the proxy an exact certificate (verify_exact.certify) is
required, otherwise the record's value is reset to the proxy.  Records already present in the
checked file (same universe, n and subset) are skipped, so the script can be rerun as results arrive."""
import json, sys, signal, os
import numpy as np
sys.path.insert(0, '/home/user/Solving-hard-math-problems-/experiments/506/search-local')
from engine import *
import run_sa as RS
from verify_exact import certify
HERE = '/home/user/Solving-hard-math-problems-/experiments/506/search-local'
OUT = f"{HERE}/runs/sa_results_checked.jsonl"
done = {}
if os.path.exists(OUT):
    for l in open(OUT):
        r = json.loads(l); done[(r['universe'], r['n'], tuple(r['S_idx']))] = r
recs = [json.loads(l) for l in open(f"{HERE}/runs/sa_results.jsonl")]
C = RS.catalogue()
cache = {}
out = open(OUT, 'a')
for r in recs:
    key = (r['universe'], r['n'], tuple(r['S_idx']))
    if key in done:
        continue
    if r['true'] >= r['proxy']:
        r['recheck'] = 'true == proxy (exact)'
        out.write(json.dumps(r) + "\n"); out.flush(); done[key] = r
        continue
    name = r['universe']
    base = name[:-2] if name.endswith('+X') else name
    NAMEMAP = {'grid4x4': 'grid4', 'grid5x5': 'grid5', 'grid6x6': 'grid6', 'grid7x7': 'grid7', 'grid4x6': 'rect4x6', 'grid3x7': 'rect3x7',
               'icosaI-vertex': 'icosa-I-v', 'icosaID-vertex': 'icosa-ID-v', 'icosaID-face': 'icosa-ID-f', 'icosaID-generic': 'icosa-ID-g',
               'icosaIDE-vertex': 'icosa-IDE-v', 'icosaIDE-generic': 'icosa-IDE-g', 'octa1-2-3-vertex': 'octa-123-v', 'octa1-2-3-generic': 'octa-123-g',
               'octa1-2-3-5-6-vertex': 'octa-12356-v', 'octa1-2-3-5-6-edge': 'octa-12356-e', 'octa1-2-3-5-6-9-vertex': 'octa-123569-v',
               'ellipse24-2-1': 'ellipse24-2-1', 'ellipse30-3-2': 'ellipse30-3-2', 'ellipse20-5-3': 'ellipse20-5-3F'}
    base = NAMEMAP.get(base, base)
    try:
        if name not in cache:
            cache[name] = RS.closure_universe(C[base](), name) if name.endswith('+X') else C[base]()
        U = cache[name]
    except Exception as e:
        r['recheck'] = f'universe rebuild failed ({e!r}); proxy value used'; r['true'] = r['proxy']; r['centre'] = None; r['deg'] = None
        out.write(json.dumps(r) + "\n"); out.flush(); done[key] = r
        print(name, r['n'], r['recheck']); continue
    S = r['S_idx']
    tc, Ot, d, nb, lines = true_count(U, S)
    note = f"recheck: old true {r['true']} -> tightened {tc}"
    if tc < r['proxy']:
        ok = False
        try:
            ex = RS.exact_coords(U, S)
            if tc > formula(r['n']) + 2:
                ex = None; note += "; far above f(n): certification skipped"
            if ex is not None:
                def _h(signum, frame): raise TimeoutError()
                signal.signal(signal.SIGALRM, _h); signal.alarm(60)
                try:
                    cert = certify(ex, Ot if isinstance(Ot, str) else [float(Ot[0]), float(Ot[1])], want=tc, verbose=False)
                finally:
                    signal.alarm(0)
                ok = cert['count'] == tc
                note += f"; exact certificate: {cert['count']}"
            else:
                note += "; no exact coordinates available"
        except Exception as e:
            note += f"; certification failed: {e!r}"[:300]
        if not ok:
            tc = r['proxy']; Ot = None; d = None
            note += "; -> proxy value used"
    print(name, r['n'], note, flush=True)
    r['true'] = int(tc); r['recheck'] = note
    r['centre'] = (None if Ot is None else (Ot if isinstance(Ot, str) else [float(Ot[0]), float(Ot[1])]))
    r['deg'] = (int(d) if d is not None else None)
    out.write(json.dumps(r) + "\n"); out.flush(); done[key] = r
out.close()
print("checked records:", len(done))
