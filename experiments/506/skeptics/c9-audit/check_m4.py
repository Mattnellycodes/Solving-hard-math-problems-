"""For enum9 output files with largest block 4 (no block of size >= 5): tabulate (D, ellmax, max degree)
and count families with no point of degree >= 8 (the only ones NOT killed by the (8_3) lemma).
Usage: python3 check_m4.py file [file ...]"""
import re, sys, collections
for fname in sys.argv[1:]:
    cnt = collections.Counter(); nodeg8 = 0; tot = 0; end = ""
    for line in open(fname):
        if line.startswith("END"):
            end = line.strip(); continue
        if not line.startswith("FAM"): continue
        tot += 1
        degs = [int(v) for v in re.search(r"degs=\[([\d,]+)\]", line).group(1).split(",")]
        D = int(re.search(r"D=(\d+)", line).group(1)); ell = int(re.search(r"ell=(\d+)", line).group(1))
        cnt[(D, ell, max(degs))] += 1
        if max(degs) < 8: nodeg8 += 1
    print(f"{fname}: {tot} labelled families with count <= 24; {end}")
    print("   (D, ellmax, maxdeg) -> count:", dict(sorted(cnt.items())))
    print("   families with NO point in >= 8 four-blocks (would survive the (8_3) lemma):", nodeg8)
