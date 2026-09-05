#!/bin/bash
# run one enumerator job with wall-clock timing: run_bg.sh n target capmode fixm tag
cd /home/user/Solving-hard-math-problems-/experiments/506/skeptics/lemmas-audit
s=$(date +%s)
./enum_own $1 $2 $3 $4 > $5.out 2> $5.err
e=$(date +%s)
echo "elapsed $((e-s)) s" >> $5.err
