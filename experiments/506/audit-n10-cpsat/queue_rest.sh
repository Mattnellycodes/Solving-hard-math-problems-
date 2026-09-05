#!/bin/bash
# remaining skeletons with largest block 5 or 6 (indices into skeletons_audit.json), sequential
for i in "$@"; do
  python3 enum_f4.py $i > f4_$i.log 2>&1
done
echo queue done "$@" >> queue_rest.done
