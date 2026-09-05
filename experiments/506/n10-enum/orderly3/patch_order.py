src = open('oe.py').read()
old = """        self.masks4 = [list_to_mask(c) for c in itertools.combinations(range(n), 4)]
        self.masks4.sort()
        self.M = len(self.masks4)
        self.pts4 = np.array([mask_to_list(m, n) for m in self.masks4], dtype=np.int8)   # (M,4)
        self.onehot"""
new = """        # codes of 4-blocks: lexicographic order of the sorted point tuples (all blocks through
        # point 0 first, then those through 1 but not 0, ...): closes the points one by one along
        # the DFS, which makes the per-point saturation bounds effective.
        self.masks4 = [list_to_mask(c) for c in itertools.combinations(range(n), 4)]
        self.M = len(self.masks4)
        self.pts4 = np.array([mask_to_list(m, n) for m in self.masks4], dtype=np.int8)   # (M,4)
        self.onehot"""
assert old in src
src = src.replace(old, new)
src = src.replace("""    Codes of 4-blocks = their index in the list of 4-subset bitmasks sorted increasingly.""",
"""    Codes of 4-blocks = their index in the lexicographically sorted list of 4-tuples.""")
open('oe.py','w').write(src)
print("patched")
