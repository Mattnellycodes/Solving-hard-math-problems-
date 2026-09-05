"""Merge search results, verify the best structures per n exactly, and write records.json.
Usage: python3 make_records.py results_*.json  [--max_alt 3]
"""
import sys, json, os, glob, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import symcore as S
import verify as V

HERE = os.path.dirname(os.path.abspath(__file__))


def describe(e, r):
    orbs = []
    k = 0; p = e['params']
    for o in e['orbits']:
        if o[0] == 'A': orbs.append(f"regular {e['m']}-gon (vertices at angles 2*pi*k/{e['m']}) radius {r['exact_params'][k]}"); k += 1
        elif o[0] == 'B': orbs.append(f"regular {e['m']}-gon (vertices at angles (2k+1)*pi/{e['m']}) radius {r['exact_params'][k]}"); k += 1
        elif o[0] == 'G': orbs.append(f"generic D_{e['m']} orbit of {2*e['m']} points, radius {r['exact_params'][k]}, angles +-({r['exact_params'][k+1]}) + 2*pi*k/{e['m']}"); k += 2
        elif o[0] == 'C': orbs.append("centre")
        elif o[0] == 'I': orbs.append("point at infinity (becomes the inversion centre after inversion)")
    s = f"D_{e['m']}-symmetric Moebius configuration: " + "; ".join(orbs) + "."
    if e['O_is_inf']:
        s += " No inversion needed (Euclidean count of the configuration itself)."
    else:
        s += f" Inverted about the point O = ({e['O'][0]:.6f}, {e['O'][1]:.6f}) (exact: intersection of the blocks {e['O_blocks'][:2]}), which lies on {e['deg']} blocks."
    s += f" Moebius blocks: {e['nblocks']} (sizes {dict((k, e['sizes'].count(k)) for k in sorted(set(e['sizes']), reverse=True))}), deg(O) = {e['deg']}."
    return s


def fmt_sizes(cs, ls):
    c = ' '.join('%d^%d' % (k, cs.count(k)) for k in sorted(set(cs), reverse=True))
    l = ' '.join('%d^%d' % (k, ls.count(k)) for k in sorted(set(ls), reverse=True)) or 'none'
    return 'circles: %s; lines: %s' % (c, l)


def main():
    files = [a for a in sys.argv[1:] if a.endswith('.json')]
    max_alt = 3
    allres = {}
    for f in files:
        d = json.load(open(f))
        for n, L in d.items():
            allres.setdefault(int(n), []).extend(L)
    records = []
    for n in sorted(allres):
        L = sorted(allres[n], key=lambda e: e['best'])
        best = L[0]['best']
        seen = set(); cnt = 0
        for e in L:
            if e['best'] > best and (e['best'] > S.formula(n) or cnt >= max_alt): continue
            if e['sig'] in seen: continue
            seen.add(e['sig'])
            print(f"== n={n} float best={e['best']} formula={S.formula(n)} {e['family']} params={e['params']}", flush=True)
            try:
                r = V.verify_entry(e, log=lambda s: print(s, flush=True))
            except Exception as ex:
                import traceback; traceback.print_exc(); r = None
            if r is None or r['degenerate']:
                print("  verification failed", flush=True); continue
            rec = {'n': n, 'circles': r['circles_exact'], 'formula_value': S.formula(n),
                   'float_best': e['best'], 'exact_matches_float': r['circles_exact'] == e['best'],
                   'construction': describe(e, r), 'family': e['family'], 'm': e['m'], 'orbits': e['orbits'],
                   'params_float': e['params'], 'params_exact': r['exact_params'],
                   'moebius_blocks': e['nblocks'], 'moebius_block_sizes': e['sizes'], 'deg_O': e['deg'],
                   'block_sizes': fmt_sizes(r['circle_sizes'], r['line_sizes']),
                   'collinear_triples': r['collinear_triples'],
                   'coordinates_exact': r['coords_sympy'], 'coordinates_float': r['coords_float'],
                   'code_path': os.path.join(HERE, 'search.py') + ' + verify.py'}
            records.append(rec); cnt += 1
    json.dump(records, open(os.path.join(HERE, 'records.json'), 'w'), indent=1)
    print("\nBEST PER n (exact):")
    bp = {}
    for rec in records:
        bp[rec['n']] = min(bp.get(rec['n'], 10**9), rec['circles'])
    for n in sorted(bp):
        print(f"n={n:2d} formula={S.formula(n):3d} best={bp[n]:3d} {'<-- BELOW' if bp[n] < S.formula(n) else ('tie' if bp[n] == S.formula(n) else 'above')}")


if __name__ == '__main__':
    main()
