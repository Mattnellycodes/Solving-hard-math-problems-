"""Lemma C (largest block) evaluation: if the largest block of an n-point Mobius configuration has m
points and r = n - m, then circles >= 1 + r*C(m,2) - m*r*(r-1)/4 - ell, ell <= (C(n,2)-o(n))/3.
Together with Lemma A (r=1 => >= f(n)) and Lemma B (r=2 => >= m^2-m+1-3*floor(m/2)) this gives, for
each n, the range of largest-block sizes m that a configuration beating f(n) must have."""
from math import comb, ceil
from mobius_enum import o_lower
def f(n): return comb(n-1,2)+1-(n-1)//2
def lemmaC(n,m):
    r=n-m
    ell=(comb(n,2)-o_lower(n))//3
    return ceil(1 + r*comb(m,2) - m*r*(r-1)/4) - ell
def lemmaB(n):
    m=n-2; return m*m-m+1-3*(m//2)
print(" n  f(n)  largest-block sizes m NOT excluded by Lemmas A,B,C (i.e. possible for count<f(n))")
for n in range(6,31):
    ok=[]
    for m in range(4,n):
        r=n-m
        if r==1: continue                      # Lemma A: count >= f(n)
        if r==2:
            if lemmaB(n)<f(n): ok.append(m)
            continue
        if lemmaC(n,m)<f(n): ok.append(m)
    print(f"{n:2d} {f(n):5d}  m in {ok}")
