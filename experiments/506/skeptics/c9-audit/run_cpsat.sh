#!/bin/bash
cd /home/user/Solving-hard-math-problems-/experiments/506/skeptics/c9-audit
python3 cpsat9_audit.py enum60 strong 1500 > cpsat_enum60_strong.log 2>&1
python3 cpsat9_audit.py m8 none 600 > cpsat_m8_none.log 2>&1
python3 cpsat9_audit.py m7 none 600 > cpsat_m7_none.log 2>&1
python3 cpsat9_audit.py m6 none 900 > cpsat_m6_none.log 2>&1
python3 cpsat9_audit.py m4_mk none 600 > cpsat_m4mk_none.log 2>&1
python3 cpsat9_audit.py m5_full none 1500 > cpsat_m5full_none.log 2>&1
python3 cpsat9_audit.py opt_full none 1800 > cpsat_optfull_none.log 2>&1
python3 cpsat9_audit.py enum60 none 1500 > cpsat_enum60_none.log 2>&1
echo DONE > run_cpsat.done
