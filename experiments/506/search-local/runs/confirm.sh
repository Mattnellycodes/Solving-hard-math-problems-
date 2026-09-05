#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-local
until grep -q "BATCH2 DONE" runs/batch2.log; do sleep 30; done
echo "===== CONFIRM START $(date)"
# confirmation rerun of the exhaustive 5x5 n=9 search with neighbour-merged arrangement clustering
python3 exhaustive_grid.py 5 9 > runs/exhaustive_grid5_n9_confirm.log 2>&1
echo "===== CONFIRM DONE $(date)"
