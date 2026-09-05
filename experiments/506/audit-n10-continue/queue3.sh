#!/bin/bash
# remaining skeletons with largest block 5 or 6 (audit numbering), one worker, sequential
cd /home/user/Solving-hard-math-problems-/experiments/506/audit-n10-continue
for i in 6 7 8 9 12 13 14 15 16 17 18 19 20 21; do
  python3 stage1_sat.py $i > runs/stage1_$i.log 2>&1
  python3 stage2_lines.py runs/stage1_$i.json > runs/stage2_$i.log 2>&1
done
python3 stage1_sat.py 11 --lb=74 --tag=_lb74 > runs/stage1_11_lb74.log 2>&1
python3 stage2_lines.py runs/stage1_11_lb74.json > runs/stage2_11_lb74.log 2>&1
python3 stage1_sat.py 11 > runs/stage1_11.log 2>&1
python3 stage2_lines.py runs/stage1_11.json > runs/stage2_11.log 2>&1
echo QUEUE3 DONE > runs/queue3.done
