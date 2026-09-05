"""Consolidate the exact certification of C-2: pass 1 (c2_exact_t*.json, norm-based sufficient
criterion) and pass 3 (c2_exact3_t*.json, embedding-aware criterion) -- and, if present, the full
pass-3 rerun over all line sets (c2_exact3_all_t*.json) as a consistency check."""
import json, os
from fractions import Fraction as Fr
def pass1_excluded(e):
    others = [iv for iv, mult in (e['delta_norm_positive_real_roots'] or []) if not (Fr(iv.strip('()').split(',')[0]) == 1 and Fr(iv.strip('()').split(',')[1]) == 1)]
    if e['reduced_gcd_degree'] == 0:
        return not others
    rr = [iv for iv, mult in e.get('reduced_gcd_norm_positive_real_roots', []) if not (Fr(iv.strip('()').split(',')[0]) == 1 and Fr(iv.strip('()').split(',')[1]) == 1)]
    return not others and not rr
status = {}
for t in ('t0', 't36'):
    for e in json.load(open(f'c2_exact_{t}.json')):
        status[(e['type_offset'], e['lineset'])] = 'EXCLUDED (pass 1: norm criterion)' if pass1_excluded(e) else 'open after pass 1'
    for e in json.load(open(f'c2_exact3_{t}.json')):
        key = (e['type_offset'], e['lineset'])
        if e['verdict'] == 'EXCLUDED':
            status[key] = status.get(key, '') if status.get(key, '').startswith('EXCLUDED') else 'EXCLUDED (pass 3: embedding-aware criterion, triple %s)' % e['triples'][-1]['triple']
n_ex = sum(1 for v in status.values() if v.startswith('EXCLUDED'))
print(f"(type, line set) cases: {len(status)}; excluded: {n_ex}; open: {[k for k, v in status.items() if not v.startswith('EXCLUDED')]}")
for k in sorted(status):
    print(f"  type {k[0]:>4} line set {k[1]:>2}: {status[k]}")
for t in ('t0', 't36'):
    fn = f'c2_exact3_all_{t}.json'
    if os.path.exists(fn):
        rep = json.load(open(fn))
        print(f"consistency rerun {fn}: {sum(1 for e in rep if e['verdict']=='EXCLUDED')} of {len(rep)} cases excluded by the embedding-aware criterion alone")
if len(status) == 62 and n_ex == 62:
    print("FINAL: no realisation of C-2 (two concentric regular pentagons, either type, any radius ratio rho != 1) has a point outside P on >= 10 of its 42 blocks; hence every planar 10-point set with block structure C-2 determines >= 42 - 9 = 33 circles.")
