#!/bin/bash
# waits for the feasibility run, then enumerates every feasible skeleton (orbit-nogood method), then post-processes
D=/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat/v2
cd $D
while pgrep -f "^python3 run_skeletons" > /dev/null; do sleep 30; done
python3 enum_skeleton.py sg 0 --time 2400 --out enum_sg_0.json > enum_sg_0.log 2>&1
for i in $(python3 -c "import json; print(' '.join(str(r['i']) for r in json.load(open('skel_sg_part2.json')) if r['status'] in ('OPTIMAL','FEASIBLE')))"); do
  python3 enum_skeleton.py sg $i --time 2400 --out enum_sg_$i.json > enum_sg_$i.log 2>&1
done
python3 postprocess.py enum_sg_*.json --tries 200 --out post_sg.json > post_sg.log 2>&1
echo QUEUE_DONE >> post_sg.log
