from itertools import combinations
n=7; pc=lambda x:bin(x).count('1')
F=[[0,1,2,3],[0,1,4,5],[2,3,4,5],[0,2,4,6],[1,3,4,6],[1,2,5,6]]
Fm=[sum(1<<i for i in b) for b in F]
triples=[sum(1<<i for i in t) for t in combinations(range(n),3)]
unc=[t for t in triples if not any(t&b==t for b in Fm)]
L=Fm+unc
sets7=[]
for S in combinations(range(len(L)),7):
    if all(pc(L[i]&L[j])<=1 for i,j in combinations(S,2)): sets7.append([L[i] for i in S])
print("line sets of size 7:",len(sets7))
for S in sets7:
    sizes=sorted(pc(x) for x in S); pairs=sum(pc(x)*(pc(x)-1)//2 for x in S)
    print("  sizes",sizes,"pairs covered",pairs,"(all 21 pairs => every pair on a >=3-point line => contradicts Sylvester-Gallai)")
