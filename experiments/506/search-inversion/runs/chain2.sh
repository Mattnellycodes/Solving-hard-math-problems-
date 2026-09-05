#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-inversion
until ! pgrep -f "[f]amily_scan.py pent" > /dev/null; do sleep 20; done
echo "===== RELAX n=11 table $(date)"
python3 relax_feas.py 11 6 table 2100 2 > runs/relax11_table.log 2>&1
echo "===== RELAX done $(date)"
for w in pent ortho polys box; do
  echo "===== SCAN $w $(date)"
  python3 family_scan.py $w > runs/family_$w.log 2>&1
done
echo "===== CHAIN DONE $(date)"
