"""Post-checks on abstract structures (F, L), own implementation.
Q = P + {inf}. Derived real point sets: at inf: P with lines L; at p in P: (P\\{p}) U {c} (c = image of inf, index n),
lines = {B\\{p} : B rich block through p, B not a line} U {(l\\{p}) U {c} : l in L, p in l}.
Tiers (applied to every subset S of every derived set, S not contained in one line):
  A: hereditary Sylvester-Gallai: S has an ordinary line (a pair of S on no restricted >=3-point line).
  B: (8_3): |S| = 8 and the restricted >=3-point lines are exactly 8 triples (Moebius-Kantor, not realisable).
  K: Kelly-Moser: ordinary lines of S >= ceil(3|S|/7).
  C: exactly-3-point restricted lines <= t3(|S|) (Burr-Gruenbaum-Sloane), ordinary lines >= o(|S|) (table)."""
import sys, json, itertools
from math import comb
from mob import popcount, bits, mask_of

T3 = {3: 1, 4: 1, 5: 2, 6: 4, 7: 6, 8: 7, 9: 10, 10: 12, 11: 16, 12: 19}
OT = {3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 4, 9: 6, 10: 5, 11: 6, 12: 6, 13: 6, 14: 7}

def derived_sets(n, blocks, lines):
    B = [mask_of(b) for b in blocks]
    L = [mask_of(l) for l in lines]
    Lset = set(L)
    out = [('inf', n, list(L))]
    c = 1 << n
    for p in range(n):
        bp = 1 << p
        ls = []
        for b in B:
            if b & bp and b not in Lset:
                ls.append(b & ~bp)
        for l in L:
            if l & bp:
                ls.append((l & ~bp) | c)
        out.append((p, n + 1, ls))
    return out

def check_set(m, lines):
    """Return dict tier -> first violation (or None). Points are 0..m-1 (some may be unused: fine)."""
    viol = {'A': None, 'B': None, 'K': None, 'C': None}
    lines = [l for l in lines if popcount(l) >= 3]
    for S in range(7, 1 << m):
        s = popcount(S)
        if s < 3:
            continue
        restr = []
        collinear = False
        for l in lines:
            r = l & S
            k = popcount(r)
            if k >= 3:
                if k == s:
                    collinear = True
                    break
                restr.append(r)
        if collinear:
            continue
        covered = sum(comb(popcount(r), 2) for r in restr)
        ordinary = comb(s, 2) - covered
        n3 = sum(1 for r in restr if popcount(r) == 3)
        if ordinary < 1 and viol['A'] is None:
            viol['A'] = (bits(S), [bits(r) for r in restr], 'no ordinary line')
        if s == 8 and len(restr) == 8 and n3 == 8 and viol['B'] is None:
            viol['B'] = (bits(S), [bits(r) for r in restr], '(8_3)')
        if ordinary < -(-3 * s // 7) and viol['K'] is None:
            viol['K'] = (bits(S), ordinary, f'ordinary {ordinary} < ceil(3*{s}/7)')
        if s in T3 and (n3 > T3[s] or ordinary < OT[s]) and viol['C'] is None:
            viol['C'] = (bits(S), n3, ordinary, f'n3 {n3} > t3({s})={T3[s]}' if n3 > T3[s] else f'ordinary {ordinary} < o({s})={OT[s]}')
    return viol

def check_structure(n, blocks, lines):
    """Return dict tier -> (where, violation) for the first violating derived set, per tier."""
    res = {'A': None, 'B': None, 'K': None, 'C': None}
    for label, m, ls in derived_sets(n, blocks, lines):
        v = check_set(m, ls)
        for t in res:
            if res[t] is None and v[t] is not None:
                res[t] = (label, v[t])
    return res

def killed_by(res, tiers):
    return [t for t in tiers if res[t] is not None]

if __name__ == '__main__':
    fn = sys.argv[1]
    d = json.load(open(fn))
    n = d['n']
    summary = []
    for i, c in enumerate(d['results']):
        blocks = c['blocks']
        rF = check_structure(n, blocks, [])
        print(f"structure {i}: sizes={[len(b) for b in blocks]} D={c['D']} ell_max={c['ell_max']} count_min={c['count_min']} degrees={c['degrees']} #line_sets={len(c['line_sets'])}")
        print(f"   blocks={blocks}")
        for t in 'ABKC':
            if rF[t] is not None:
                print(f"   F alone, tier {t}: violated at {rF[t][0]}: {rF[t][1]}")
        surv = {'A': [], 'AB': [], 'ABK': [], 'ABKC': []}
        for j, L in enumerate(c['line_sets']):
            r = check_structure(n, blocks, L)
            cnt = comb(n, 3) - c['D'] - len(L)
            if r['A'] is None:
                surv['A'].append((j, cnt))
                if r['B'] is None:
                    surv['AB'].append((j, cnt))
                    if r['K'] is None:
                        surv['ABK'].append((j, cnt))
                        if r['C'] is None:
                            surv['ABKC'].append((j, cnt))
        for k, v in surv.items():
            print(f"   line sets surviving tier {k}: {len(v)} counts={sorted(set(x[1] for x in v))}")
        for (j, cnt) in surv['ABK']:
            print(f"      A+B+K survivor: line set {j} count={cnt} lines={c['line_sets'][j]}  passes C: {(j,cnt) in surv['ABKC']}")
        summary.append((i, [len(b) for b in blocks].count(5), {k: len(v) for k, v in surv.items()}))
    print("SUMMARY", summary)
