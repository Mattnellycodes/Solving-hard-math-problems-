"""Locate the two example structures reported by n10-cpsat/v2 (skeleton 0 count-31, skeleton 9 count-32)
in the auditor's (F,L) class lists, and report which filters kill them."""
import json
from common import incidence_graph, isomorphic, run_filters, count_of, check_structure
v2 = {
 'v2 skeleton 0 example (6-block, count 31)': (
   [[0,1,2,3,4,5],[0,1,6,8],[0,1,7,9],[0,2,6,9],[0,3,8,9],[0,4,7,8],[0,5,6,7],[1,2,8,9],[1,3,6,7],[1,4,6,9],[1,5,7,8],[2,3,7,8],[2,4,6,7],[2,5,6,8],[2,5,7,9],[3,4,6,8],[3,4,7,9],[3,5,6,9],[4,5,8,9],[6,7,8,9]],
   [[0,2,7],[0,3,6],[0,4,9],[0,5,8],[1,2,6],[1,3,8],[1,4,7],[1,5,9],[2,3,9],[2,4,8],[3,5,7],[4,5,6],[6,7,8,9]], 'fl_1.json'),
 'v2 skeleton 9 example (two disjoint 5-blocks, count 32)': (
   [[0,1,2,3,4],[0,1,5,6],[0,1,7,9],[0,2,5,8],[0,2,6,7],[0,3,5,7],[0,3,8,9],[0,4,5,9],[0,4,6,8],[1,2,5,9],[1,2,6,8],[1,3,5,8],[1,3,6,7],[1,4,6,9],[1,4,7,8],[2,3,6,9],[2,3,7,8],[2,4,5,7],[2,4,8,9],[3,4,5,6],[3,4,7,9],[5,6,7,8,9]],
   [[0,2,9],[0,3,6],[0,7,8],[1,2,7],[1,4,5],[1,8,9],[2,5,6],[3,4,8],[3,5,9],[4,6,7]], 'fl_7.json'),
}
for name, (F, L, src) in v2.items():
    check_structure(F, L)
    G = incidence_graph(F, L)
    cls = json.load(open(src))['classes']
    hits = [i for i, c in enumerate(cls) if isomorphic(G, incidence_graph(c['F'], c['L']))]
    print(name, 'count', count_of(F, L), '-> isomorphic to class(es)', hits, 'of', src)
    print('   killed by:', sorted(set(k for _, k, _ in run_filters(F, L))))
