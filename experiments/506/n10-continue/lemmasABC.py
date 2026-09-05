"""Re-derived bounds excluding a largest block of size m >= 7 for n = 10 (target: circles <= 32).
A (m = 9): blocks through x: the 36 triples {x,a,b} lie in distinct blocks; l <= 4 (B a circle: lines through x cut B
   in disjoint pairs) or l = 1 (B a line: no other line): circles >= 37 - 4 = 33 = f(10).
B (m = 8): |blocks| = 65 - 3t with t <= 4 four-blocks {x,y,a,b}: >= 53; l <= 14 (SG: lines cover <= 44 pairs) => >= 39.
C (m = 7): chunk lemma N2 >= r*C(m,2) - m*r(r-1)/4 = 63 - 10.5 -> >= 53 blocks meeting B in 2 points; +B itself;
   l <= 14 => circles >= 54 - 14 = 40.
All three exceed 32, so a configuration with <= 32 circles has largest block of size <= 6 (and >= 4: without rich
blocks D = 0 and l <= 14 gives circles >= 106)."""
import math
n = 10
for m in (9, 8, 7, 6, 5):
    r = n - m
    lmax = (math.comb(n, 2) - 1) // 3
    if m == 9:
        print("m=9 (Lemma A): circles >= 1 + 36 - 4 =", 37 - 4)
    elif m == 8:
        print("m=8 (Lemma B): blocks >= 65 - 12 = 53, circles >= 53 -", lmax, "=", 53 - lmax)
    else:
        N2 = math.ceil(r * math.comb(m, 2) - m * r * (r - 1) / 4)
        print(f"m={m} (Lemma C): N2 >= {N2}, circles >= {1 + N2} - {lmax} = {1 + N2 - lmax}")
