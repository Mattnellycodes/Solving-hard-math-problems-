#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/search-local
nohup bash runs/batchA.sh > runs/batchA.log 2>&1 &
nohup bash runs/batchB.sh > runs/batchB.log 2>&1 &
