"""Isomorphism-class dedupe of FAM lines produced by enum9 (from scratch).
Canonical form of a block family on n points = lexicographically smallest sorted tuple of block
bitmasks over ALL n! relabellings (vectorised with numpy; n! = 362880 for n = 9).
Usage: python3 iso.py file1.txt [file2.txt ...]
"""
import sys, re, itertools, ast
import numpy as np
from math import comb

def parse(fname):
    fams = []
    for line in open(fname):
        if line.startswith("FAM"):
            mD = re.search(r"D=(\d+)", line); ml = re.search(r"ell=(\d+)", line); mc = re.search(r"count=(-?\d+)", line)
            blocks = ast.literal_eval(line.split("blocks=")[1].strip())
            fams.append(dict(D=int(mD.group(1)), ell=int(ml.group(1)), count=int(mc.group(1)), blocks=[tuple(b) for b in blocks]))
        elif line.startswith("END"):
            fams.append(dict(END=line.strip()))
    return fams

_perm_cache = {}
def perm_table(n):
    if n not in _perm_cache:
        perms = np.array(list(itertools.permutations(range(n))), dtype=np.int64)   # (n!, n)
        _perm_cache[n] = (1 << perms)                                               # images of singletons
    return _perm_cache[n]

def canon(blocks, n):
    """blocks: list of tuples of point labels.  Returns canonical sorted tuple of codes."""
    Pw = perm_table(n)                                   # (n!, n): Pw[s, i] = 2^{s(i)}
    imgs = np.zeros((Pw.shape[0], len(blocks)), dtype=np.int64)
    for j, b in enumerate(blocks):
        imgs[:, j] = Pw[:, list(b)].sum(axis=1)
    imgs.sort(axis=1)
    # lexicographic minimum over rows
    order = np.lexsort(imgs.T[::-1])
    return tuple(int(v) for v in imgs[order[0]])

def decode(code, n):
    return tuple(i for i in range(n) if code >> i & 1)

def main():
    n = 9
    for fname in sys.argv[1:]:
        fams = parse(fname)
        ends = [f["END"] for f in fams if "END" in f]
        fams = [f for f in fams if "END" not in f]
        classes = {}
        for f in fams:
            c = canon(f["blocks"], n)
            classes.setdefault(c, []).append(f)
        print(f"== {fname}: {len(fams)} labelled families, {len(classes)} isomorphism classes; {ends}")
        for c, lst in classes.items():
            f = lst[0]
            blocks = [decode(x, n) for x in c]
            degs = [sum(1 for b in blocks if p in b) for p in range(n)]
            sizes = sorted((len(b) for b in blocks), reverse=True)
            print(f"   class: {len(lst)} labelled copies; sizes={sizes} D={f['D']} ellmax={f['ell']} count={f['count']} degrees={degs}")
            print(f"      canonical blocks: {[list(b) for b in blocks]}")

if __name__ == "__main__":
    main()
