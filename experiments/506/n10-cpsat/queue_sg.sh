#!/bin/bash
D=/home/user/Solving-hard-math-problems-/experiments/506/n10-cpsat
while pgrep -f "run_case.py A s[g]" > /dev/null; do sleep 10; done
cd $D && python3 extend.py sg --orchard10 --time 1500 --out extend_sg.json > extend_sg.log 2>&1
