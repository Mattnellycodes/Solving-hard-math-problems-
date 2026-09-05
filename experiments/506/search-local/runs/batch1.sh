#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-local
python3 run_sa.py grid4 grid5 grid6 rect4x6 rect3x7 grid5inv-c grid5inv-corner grid4inv-multi grid5inv-multi grid6inv-c trilat3 trilat4 trilat5 hyperbola14 hyperbola23 --ns 9-16 --iters 300000 --restarts 6 --seed 1
python3 run_sa.py poly5-phi poly8-s2 poly6-s3 poly12 poly7 poly9 poly10 poly14 poly16 poly18 poly20 poly24 poly30 --ns 9-16 --iters 300000 --restarts 6 --seed 2
echo BATCH1 DONE
