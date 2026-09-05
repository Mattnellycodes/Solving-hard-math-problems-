#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/skeptics/c9-audit
until grep -q "done m=4" enum9_none_m4.err 2>/dev/null || ! pgrep -f "enum9 4 28 36" >/dev/null; do sleep 10; done
python3 check_m4.py enum9_strong_m4.txt enum9_sg_m4.txt enum9_none_m4.txt > check_m4.log 2>&1
python3 angle_lemma.py > angle_lemma.log 2>&1
echo AFTER_M4_DONE >> check_m4.log
