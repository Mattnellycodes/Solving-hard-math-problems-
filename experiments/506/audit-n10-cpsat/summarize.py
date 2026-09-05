"""Overview table of the audit results (reads f4_*.json, fl_*.json, surv_*.json)."""
import json, os, glob
sk = json.load(open('skeletons_audit.json'))
print('| idx | big blocks | status | labelled | F-classes | (F,L) classes | survivors of proved filters |')
print('|---|---|---|---|---|---|---|')
for i, F in enumerate(sk):
    if any(len(B) >= 7 for B in F):
        print('| %d | %s | excluded (block >= 7, Lemmas A/B/C) | | | | |' % (i, sorted(len(B) for B in F))); continue
    f4 = 'f4_%d.json' % i
    if not os.path.exists(f4):
        print('| %d | %s | (not run) | | | | |' % (i, sorted(len(B) for B in F))); continue
    d = json.load(open(f4))
    nf = len(d['classes'])
    fl = 'fl_%d.json' % i
    nfl = len(json.load(open(fl))['classes']) if os.path.exists(fl) else (0 if nf == 0 else '?')
    sv = 'surv_%d.json' % i
    ns = len(json.load(open(sv))) if os.path.exists(sv) else (0 if nfl == 0 else '?')
    print('| %d | %s | %s | %d | %d | %s | %s |' % (i, sorted(len(B) for B in F), d['status'], d['labelled'], nf, nfl, ns))
if os.path.exists('f4_empty.json'):
    d = json.load(open('f4_empty.json'))
    print('| - | [] (all blocks of size <= 4) | %s | %d | %d | | |' % (d['status'], d['labelled'], len(d['classes'])))
