#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/n10-continue
for i in 7 5 6 8 9 11 12 13 14 15 16 17 18 19 20 21; do
  python3 enum4.py $i --out runs/enum_$i.json > runs/enum_$i.log 2>&1
done
echo ALL_DONE > runs/rest_done.flag
