"""No-lex validation for the single-6-block skeleton: the lex-leader run found 8 F-classes; the run
without symmetry breaking found 9840 labelled solutions.  If sum_i |Aut(skel)| / |Aut(F_i)| over the
8 classes equals 9840, the no-lex solution set is exactly the union of the 8 orbits (each orbit is
contained in it), so no class is missing."""
import json
from common import automorphisms
d = json.load(open('f4_1.json'))
skel = [B for B in d['skeleton']]
A = len(automorphisms(skel, []))
tot = 0
for cl in d['classes']:
    a = len(automorphisms(cl['F'], []))
    tot += A // a
    print('|Aut(F)| =', a, 'orbit size', A // a)
print('|Aut(skeleton)| =', A, ' sum of orbit sizes =', tot, ' no-lex labelled solutions = 9840 ->', 'MATCH' if tot == 9840 else 'MISMATCH')
