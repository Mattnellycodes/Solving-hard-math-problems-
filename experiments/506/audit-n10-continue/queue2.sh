#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/audit-n10-continue
python3 stage1_sat.py 0 --nocap --tag=_nocap > runs/stage1_0_nocap.log 2>&1
python3 stage1_sat.py 0 --lb=74 --tag=_lb74 > runs/stage1_0_lb74.log 2>&1
python3 stage1_sat.py 0 --solver=glucose4 --tag=_glucose > runs/stage1_0_glucose.log 2>&1
python3 stage1_sat.py 1 --nocap --tag=_nocap > runs/stage1_1_nocap.log 2>&1
python3 stage1_sat.py 5 --nocap --tag=_nocap > runs/stage1_5_nocap.log 2>&1
echo DONE > runs/queue2.done
