#!/bin/bash
# runs after batch1.sh has finished: the conic universes
cd /home/user/Solving-hard-math-problems-/experiments/506/search-inversion
until ! pgrep -f "^bash runs/batch1.sh" > /dev/null; do sleep 30; done
echo "===== GROUP conic $(date)"
python3 run_all.py conic 9 16 10 3 12
echo "===== BATCH2 DONE $(date)"
