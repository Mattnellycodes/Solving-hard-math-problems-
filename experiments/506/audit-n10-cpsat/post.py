"""Apply the auditor's filters (common.run_filters) to every (F,L) class of an fl_<tag>.json."""
import sys, json, time
from common import run_filters, count_of
for src in sys.argv[1:]:
    data = json.load(open(src))
    t0 = time.time()
    print('==', src, len(data['classes']), 'classes')
    surv = []
    for i, cl in enumerate(data['classes']):
        F, L = cl['F'], cl['L']
        v = run_filters(F, L)
        vc = run_filters(F, L, cited=True)
        kinds = sorted(set(k for _, k, _ in v))
        kindsc = sorted(set(k for _, k, _ in vc))
        print('class %d count=%d l=%d: proved-filters -> %s ; +cited -> %s' % (i, count_of(F, L), len(L), kinds or 'SURVIVES', kindsc or 'SURVIVES'))
        for name, k, w in v[:3]:
            print('     e.g.', name, k, w)
        if not v:
            surv.append(cl)
    print('survivors (proved filters):', len(surv), ' (%.0fs)' % (time.time() - t0))
    json.dump(surv, open('surv_' + src[3:], 'w'))
