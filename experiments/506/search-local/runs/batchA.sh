#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-local
echo "===== BATCHA START $(date)"
python3 run_sa.py grid5inv-corner grid4inv-multi grid5inv-multi grid6inv-c trilat3 trilat4 trilat5 hyperbola14 hyperbola23 --ns 9-16 --iters 300000 --restarts 6 --seed 1 --skip-done
python3 run_sa.py poly5-phi poly8-s2 poly6-s3 poly12 poly7 poly9 poly10 poly14 poly16 poly18 poly20 poly24 poly30 --ns 9-16 --iters 300000 --restarts 6 --seed 2 --skip-done
python3 run_sa.py grid4 ortho-A1 icosa-I-v octa-123-v poly8-s2 poly6-s3 poly5-phi hyperbola14 ellipse24-2-1 --ns 9-16 --iters 200000 --restarts 4 --seed 7 --closure --skip-done
echo "===== BATCHA DONE $(date)"
