from math import comb
def f(n): return comb(n-1,2)+1-(n-1)//2
t3 = {3:1,4:1,5:2,6:4,7:6,8:7,9:10,10:12,11:16,12:19}   # Burr-Grunbaum-Sloane 1974 (quoted by the report)
print("n, f(n), needed D+l >= C(n,3)-f(n)+1, claimed:")
claimed = dict(zip(range(9,17),[60,88,125,170,226,292,371,462]))
for n in range(9,17):
    need = comb(n,3) - f(n) + 1
    print(n, f(n), need, "== C(n-1,3)+floor((n-1)/2) =", comb(n-1,3)+(n-1)//2, "claimed", claimed[n], "OK" if need==claimed[n]==comb(n-1,3)+(n-1)//2 else "BAD")
print("all-4-block cap 3*n*t3(n-1)/4 + t3(n):")
for n in range(9,13):
    cap = (3*n*t3[n-1])//4 + t3[n]
    print(n, cap, "claimed", {9:57,10:87,11:115,12:163}[n], "need", comb(n-1,3)+(n-1)//2)
