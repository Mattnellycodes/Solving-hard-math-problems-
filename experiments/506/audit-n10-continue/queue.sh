#!/bin/bash
# remaining skeletons, one worker, sequential
cd /home/user/Solving-hard-math-problems-/experiments/506/audit-n10-continue
for i in 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 2 3 10; do
  python3 stage1_sat.py $i > runs/stage1_$i.log 2>&1
  python3 stage2_lines.py runs/stage1_$i.json > runs/stage2_$i.log 2>&1
done
echo QUEUE DONE > runs/queue.done
