import json, itertools
from mob import *
d = json.load(open('out_n10_t32_sg.json')); PT = PermTable(10)
for i in (4, 5, 6, 7, 8):
    masks = sorted(mask_of(b) for b in d['results'][i]['blocks'])
    less, equal = rows_less_or_equal(PT.image_masks(masks), masks)
    print(f"structure {i}: |Aut(F)| = {int(equal.sum())}")
# (8_3): every family of 8 triples on 8 points pairwise sharing <= 1 point is one isomorphism class
PT8 = PermTable(8)
trip = [mask_of(t) for t in itertools.combinations(range(8), 3)]
fams = []
def rec(F, start):
    if len(F) == 8: fams.append(list(F)); return
    for i in range(start, len(trip)):
        if all(popcount(trip[i] & b) <= 1 for b in F): rec(F + [trip[i]], i + 1)
rec([], 0)
cl = set(canon_full(PT8, [bits(m) for m in F]) for F in fams)
print(f"(8_3): {len(fams)} labelled families of 8 pairwise <=1-intersecting triples on 8 points, {len(cl)} isomorphism class(es); 8!/48 = {40320//48}")
