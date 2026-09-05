#!/bin/sh
# Reproduce everything (total ~7 minutes on one core; the exact k=7 brute force dominates).
set -e
cd "$(dirname "$0")"
mkdir -p logs
python3 part1_exact_small_n.py      | tee logs/part1.txt
python3 route_a_iv.py               | tee logs/route_a.txt
python3 route_b_exact.py            | tee logs/route_b.txt
python3 check_certificate.py cert_route_a.json cert_route_b.json | tee logs/checker.txt
python3 crosscheck_bruteforce.py    | tee logs/crosscheck_bruteforce.txt
python3 validate_lemmas.py          | tee logs/validate_lemmas.txt
python3 uniform_bound.py            | tee logs/uniform_bound.txt
