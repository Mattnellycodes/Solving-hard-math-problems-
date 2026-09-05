#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-inversion
# wait for the pentagon scan to finish
until ! pgrep -f "family_scan.py pent" > /dev/null; do sleep 20; done
for w in ortho polys box; do
  echo "===== SCAN $w $(date)"
  python3 family_scan.py $w
done
echo "===== SCANS DONE $(date)"
