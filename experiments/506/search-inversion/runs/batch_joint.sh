#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-inversion
export JOINT=1
for g in lattice poly tri polygon grid conic pencil klein; do
  echo "===== GROUP $g $(date)"
  python3 run_all.py $g 9 16 60 0 12
done
echo "===== JOINT BATCH DONE $(date)"
