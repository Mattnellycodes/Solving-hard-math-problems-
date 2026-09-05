#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/skeptics/c6-construct
: > runner.log
for f in grid10 tri; do
  echo "== $f (exact)" >> runner.log
  ./brute6 < $f.txt > ${f}_out.txt 2>> runner.log
done
for p in pool_*.txt; do
  n=${p%.txt}
  echo "== $n (float)" >> runner.log
  ./brute6f < $p > ${n}_out.txt 2>> runner.log
done
echo "ALL DONE" >> runner.log
