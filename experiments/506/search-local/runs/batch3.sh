#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-local
until grep -q "BATCH1 DONE" runs/batch1.log; do sleep 30; done
echo "===== BATCH3 START $(date)"
# level-2 orthocentric closures (exact block finder)
python3 run_sa.py ortho-A2 ortho-B2 ortho-C2 ortho-D2 ortho-E2 ortho-F1 ortho-G2 ortho-H2 --ns 9-16 --iters 200000 --restarts 4 --seed 3
# coincidence-closed universes (+X)
python3 run_sa.py grid4 ortho-A1 icosa-I-v octa-123-v poly8-s2 poly6-s3 poly5-phi hyperbola14 ellipse24-2-1 --ns 9-16 --iters 200000 --restarts 4 --seed 7 --closure
# second seeds on the richest universes
python3 run_sa.py grid5inv-c grid5inv-multi trilat4 poly12 poly10 --ns 9-16 --iters 400000 --restarts 6 --seed 8
echo "===== BATCH3 DONE $(date)"
