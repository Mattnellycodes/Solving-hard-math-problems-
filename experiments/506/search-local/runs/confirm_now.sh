#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-local
B=$(ps -eo pid,args | grep "python3 run_s[a].py ortho-A2" | awk '{print $1}')
[ -n "$B" ] && kill -STOP $B
echo "paused batch B pid=$B; CONFIRM START $(date)"
python3 exhaustive_grid.py 5 9 > runs/exhaustive_grid5_n9_confirm.log 2>&1
echo "CONFIRM DONE $(date)"
[ -n "$B" ] && kill -CONT $B
echo "resumed batch B"
