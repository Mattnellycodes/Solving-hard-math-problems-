"""(V2) Independent verification of the n=8 case analysis WITHOUT the orderly-generation machinery.

Facts used (proved in REPORT.md):
  * count = 56 - D - ell;  count <= 17 forces max block size 4 and b4 in {11,12}
    (5-,6-,7-block cases are killed by counting; see REPORT Prop. 4.1) -- here we do NOT assume that
    and instead check all families with a 5-,6- or 7-block by brute force too.
  * every point lies in <= 6 four-blocks (7 would give a derived Fano plane: 0 ordinary lines).
  * b4 = 12 => all degrees 6;  b4 = 11 => >= 4 points of degree 6.
  * 6 four-blocks through p with pairwise <=1 common point in P\{p} form 'Fano minus a line'
    (unique up to isomorphism), so WLOG the blocks through point 0 are
    {0}+{1,2,4},{0}+{2,3,5},{0}+{3,4,6},{0}+{4,5,7},{0}+{5,6,1},{0}+{6,7,2}   (Fano lines {i,i+1,i+3} mod 7
    minus the line {7,1,3}).
This script: enumerates all 4-block families containing those 6 blocks with 11 or 12 blocks, checks
degrees, computes ell_max by brute force (all subsets of pairwise <=1-intersecting blocks among
4-blocks and uncovered triples, with pair budget <= 24), and reports min count.
Also brute-forces the cases with a block of size >= 5 by direct labelled DFS (small).
"""
import itertools
from math import comb
n=8
def pc(x): return bin(x).count('1')
def code(s): return sum(1<<i for i in s)
fano=[{1,2,4},{2,3,5},{3,4,6},{4,5,7},{5,6,1},{6,7,2},{7,1,3}]
base=[code({0}|l) for l in fano[:6]]
all4=[c for c in range(1<<n) if pc(c)==4]
def compatible(F,b): return all(pc(a&b)<=2 for a in F)
def ell_max(F):
    unc=[code(t) for t in itertools.combinations(range(n),3) if not any(code(t)&b==code(t) for b in F)]
    cands=sorted(F)+unc
    best=0
    def rec(i,chosen,pairs):
        nonlocal best
        best=max(best,len(chosen))
        for j in range(i,len(cands)):
            c=cands[j]
            if pairs+comb(pc(c),2)>24: continue
            if all(pc(c&d)<=1 for d in chosen): rec(j+1,chosen+[c],pairs+comb(pc(c),2))
    rec(0,[],0)
    return best
results={}
def rec(F,start):
    k=len(F)
    if k in (11,12):
        deg=[sum(1 for b in F if b>>p&1) for p in range(n)]
        if max(deg)<=6:
            e=ell_max(F); D=3*k; cnt=56-D-e
            results.setdefault((k,e),[]).append(sorted(F))
    if k==12: return
    for j in range(start,len(all4)):
        b=all4[j]
        if b in F or not compatible(F,b): continue
        if b>>0&1: continue          # point 0 already has degree 6
        deg_ok=all(sum(1 for a in F if a>>p&1)<6 for p in range(n) if b>>p&1)
        if not deg_ok: continue
        rec(F+[b],j+1)
rec(base,0)
for (k,e),fams in sorted(results.items()):
    print(f"b4={k}: ell_max={e}: {len(fams)} labelled families, count = {56-3*k-e}")
# classify the 12-block families up to isomorphism (brute force canonical form)
perms=list(itertools.permutations(range(n)))
def img(c,s): return sum(1<<s[i] for i in range(n) if c>>i&1)
def canon(F): return min(tuple(sorted(img(c,s) for c in F)) for s in perms)
twelve=[F for (k,e),fams in results.items() if k==12 for F in fams]
print("12-block families (labelled, containing base):",len(twelve),"iso classes:",len(set(canon(F) for F in twelve)))
# check they are all AG(3,2) minus a parallel class: i.e. adding two complementary 4-sets gives an SQS(8)
def is_sqs(F):
    return len(F)==14 and all(pc(a&b)<=2 for a,b in itertools.combinations(F,2))
ok=0
for F in twelve:
    found=False
    for a in all4:
        b=(~a)&255
        if a<b and a not in F and b not in F and compatible(F,a) and compatible(F,b) and is_sqs(F+[a,b]): found=True
    ok+=found
print("12-block families completable to SQS(8) by a parallel class:",ok,"of",len(twelve))
# Brute force: families containing a block of size >=5 (any), count<=17 possible?  D + ell >= 39.
big=[c for c in range(1<<n) if 5<=pc(c)<=7]
best=None
def rec2(F,last,D):
    global best
    if D>=31:   # ell<=8
        e=ell_max(F)
        if 56-D-e<=17: print("!! big-block candidate",[sorted(i for i in range(n) if b>>i&1) for b in F],56-D-e)
    cands=[c for c in range(1<<n) if 4<=pc(c)<=7 and c>last and compatible(F,c)]
    pot=sum(comb(pc(c),3)-1 for c in cands)
    if D+pot<31: return
    for c in cands:
        # coverage cap: each point covered <= 18 pairs by big blocks
        cov=[sum(comb(pc(a)-1,2) for a in F+[c] if a>>p&1) for p in range(n)]
        if max(cov)>18: continue
        rec2(F+[c],c,D+comb(pc(c),3)-1)
# start with one big block (by symmetry {0..k-1}) then extend with anything
for k in (5,6,7):
    b=code(range(k)); rec2([b],-1,comb(k,3)-1)
print("big-block brute force done (no output above => none reach 17)")
