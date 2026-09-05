#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-inversion
for g in lattice poly tri polygon grid pencil klein; do
  echo "===== GROUP $g $(date)"
  python3 run_all.py $g 9 16 10 3 12
done
echo "===== BATCH DONE $(date)"
