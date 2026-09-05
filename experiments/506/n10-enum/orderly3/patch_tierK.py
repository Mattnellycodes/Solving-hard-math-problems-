src = open('hcheck.py').read()
old = """        if 'C' in tiers:
            n3 = sum(1 for x in RS if popcount(x) == 3)"""
new = """        if 'K' in tiers and s >= 3 and ordinary < -(-3 * s // 7):     # Kelly-Moser: o(m) >= 3m/7
            viol['K'].append(f"subset {mask_to_list(S, m)} has {ordinary} ordinary lines < ceil(3*{s}/7)={-(-3 * s // 7)}")
        if 'C' in tiers:
            n3 = sum(1 for x in RS if popcount(x) == 3)"""
assert old in src; src = src.replace(old, new)
# default tiers everywhere: A, B, K, C
src = src.replace("tiers=('A', 'B', 'C')", "tiers=('A', 'B', 'K', 'C')")
old2 = """        for t in ('A', 'B', 'C'):
            fv = rep['F_only'][t]"""
new2 = """        for t in ('A', 'B', 'K', 'C'):
            fv = rep['F_only'][t]"""
assert old2 in src; src = src.replace(old2, new2)
old3 = """        print(f"   line sets surviving: tier A: {len(cum['A'])} (counts {counts['A']}); "
              f"A+B: {len(cum['B'])} (counts {counts['B']}); A+B+C: {len(cum['C'])} (counts {counts['C']})")
        if rep['line_sets']:
            # show one example violation for the first line set for each tier
            lv = rep['line_sets'][0]
            for t in ('A', 'B', 'C'):"""
new3 = """        print(f"   line sets surviving: tier A: {len(cum['A'])} (counts {counts['A']}); "
              f"A+B: {len(cum['B'])} (counts {counts['B']}); A+B+K: {len(cum['K'])} (counts {counts['K']}); "
              f"A+B+K+C: {len(cum['C'])} (counts {counts['C']})")
        if rep['line_sets']:
            # show one example violation for the first line set for each tier
            lv = rep['line_sets'][0]
            for t in ('A', 'B', 'K', 'C'):"""
assert old3 in src; src = src.replace(old3, new3)
old4 = """        verdict = ("SURVIVES all tiers" if cum['C'] else
                   "killed by tier C only (cited tables)" if cum['B'] else"""
new4 = """        verdict = ("SURVIVES all tiers" if cum['C'] else
                   "killed by tier C only (orchard / exact o-table)" if cum['K'] else
                   "killed by tier K (Kelly-Moser o(m) >= 3m/7)" if cum['B'] else"""
assert old4 in src; src = src.replace(old4, new4)
old5 = """        summary.append((i, r['sizes'].count(5), r['sizes'].count(6), len(r['blocks']), r['count_min'],
                        len(cum['A']), len(cum['B']), len(cum['C']), counts, verdict))
    print("\\nSUMMARY (structure, #5-blocks, #6-blocks, #blocks, count_min, surviving line sets A / A+B / A+B+C, counts, verdict)")"""
new5 = """        summary.append((i, r['sizes'].count(5), r['sizes'].count(6), len(r['blocks']), r['count_min'],
                        len(cum['A']), len(cum['B']), len(cum['K']), len(cum['C']), counts, verdict))
    print("\\nSUMMARY (structure, #5-blocks, #6-blocks, #blocks, count_min, surviving line sets A / A+B / A+B+K / A+B+K+C, counts, verdict)")"""
assert old5 in src; src = src.replace(old5, new5)
old6 = """        json.dump([{"index": s[0], "verdict": s[9], "surv_A": s[5], "surv_AB": s[6], "surv_ABC": s[7],
                    "counts": {k: v for k, v in s[8].items()}} for s in summary], open(sys.argv[2], 'w'), indent=1)"""
new6 = """        json.dump([{"index": s[0], "verdict": s[10], "surv_A": s[5], "surv_AB": s[6], "surv_ABK": s[7],
                    "surv_ABKC": s[8], "counts": {k: v for k, v in s[9].items()}} for s in summary],
                  open(sys.argv[2], 'w'), indent=1)"""
assert old6 in src; src = src.replace(old6, new6)
open('hcheck.py', 'w').write(src)
print("patched")
