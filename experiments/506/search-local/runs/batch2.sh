#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-local
# wait for the exhaustive job to finish
# (wait loop removed)
echo "===== BATCH2 START $(date)"
python3 run_sa.py ortho-A1 ortho-A2 ortho-B2 ortho-C2 ortho-D2 ortho-E2 ortho-F1 ortho-G2 ortho-H2 --ns 9-16 --iters 300000 --restarts 6 --seed 3
python3 run_sa.py icosa-I-v icosa-ID-v icosa-ID-f icosa-ID-g icosa-IDE-v icosa-IDE-g octa-123-v octa-123-g octa-12356-v octa-12356-e octa-123569-v --ns 9-16 --iters 300000 --restarts 6 --seed 4
python3 run_sa.py ellipse24-2-1 ellipse24-2-1F ellipse30-3-2 ellipse20-5-3F ellipse24-2-1-rot hyperbola-pow --ns 9-16 --iters 300000 --restarts 6 --seed 5
python3 run_sa.py pdiag8 pdiag10 pdiag12 pdiag16 pdiag18 pdiag24 pdiag30 pdiag8x2 pdiag12x2 --ns 9-16 --iters 150000 --restarts 4 --seed 6
echo "===== BATCH2 DONE $(date)"
