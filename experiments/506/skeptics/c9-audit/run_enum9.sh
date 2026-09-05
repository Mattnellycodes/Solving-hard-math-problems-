#!/bin/bash
# n=9, target 24, three cap regimes; largest block m = 5..8 (complete), then m = 4 with a node limit.
cd /home/user/Solving-hard-math-problems-/experiments/506/skeptics/c9-audit
for regime in "strong 24 30" "sg 27 35" "none 28 36"; do
  set -- $regime; name=$1; cp=$2; cl=$3
  for m in 5 6 7 8; do
    ./enum9 $m $cp $cl 24 > enum9_${name}_m${m}.txt 2> enum9_${name}_m${m}.err
  done
done
for regime in "strong 24 30" "sg 27 35" "none 28 36"; do
  set -- $regime; name=$1; cp=$2; cl=$3
  ./enum9 4 $cp $cl 24 3000000000 > enum9_${name}_m4.txt 2> enum9_${name}_m4.err
done
echo ALLDONE > run_enum9.done
