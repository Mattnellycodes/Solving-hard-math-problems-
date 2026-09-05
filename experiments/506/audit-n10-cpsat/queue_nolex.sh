#!/bin/bash
python3 enum_f4.py 0 --nolex > f4_0_nolex.log 2>&1
python3 enum_f4.py 1 --nolex > f4_1_nolex.log 2>&1
echo nolex done >> queue_rest.done
