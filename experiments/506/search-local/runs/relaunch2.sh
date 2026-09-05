#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-local
python3 - <<'PY'
import json
HERE='/home/user/Solving-hard-math-problems-/experiments/506/search-local'
keep, arch = [], []
for line in open(f"{HERE}/runs/sa_results.jsonl"):
    try: r = json.loads(line)
    except Exception: continue
    (arch if (r['universe'] == 'grid4inv-multi' and r['n'] >= 15) or r['universe'] == 'ortho-F1' else keep).append(line)
open(f"{HERE}/runs/sa_results.jsonl", 'w').writelines(keep)
open(f"{HERE}/runs/sa_results_archived_lowbudget.jsonl", 'a').writelines(arch)
print("kept", len(keep), "archived", len(arch))
PY
# drop the (already done) exhaustive confirmation from batch B
sed -i '/exhaustive_grid.py 5 9/d; /CONFIRM START/d; /CONFIRM DONE/d; /confirmation rerun/d' runs/batchB.sh
echo "===== RELAUNCH $(date)" >> runs/batchA.log; echo "===== RELAUNCH $(date)" >> runs/batchB.log
nohup bash runs/batchA.sh >> runs/batchA.log 2>&1 &
nohup bash runs/batchB.sh >> runs/batchB.log 2>&1 &
sleep 3
ps -eo pid,args | grep "run_s[a]" | cut -c1-90
