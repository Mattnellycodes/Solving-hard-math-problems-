# Part 3: (a) Binet integral representation + remainder factor theta in (0,1); (b) non-squarefree formula vs exact rationals;
# (c) refine the q=2 threshold; (d) directed-rounding certificate for UB(92) from exact rationals.
from fractions import Fraction
from mpmath import mp, mpf, euler, pi, quad, exp, log, psi, bernoulli, fprod, harmonic, iv, zeta
from sympy import primerange, divisors, mobius, factorint
mp.dps = 40
# (a) psi(z) = ln z - 1/(2z) - 2 int_0^inf t/((t^2+z^2)(e^{2 pi t}-1)) dt   (Binet, differentiated)
for z in [mpf(1), mpf(2), mpf('3.5'), mpf(30)]:
    I = quad(lambda t: t/((t*t+z*z)*(exp(2*pi*t)-1)), [0, mp.inf])
    print("z=",z," psi(z)=",psi(0,z)," Binet rhs=", log(z)-1/(2*z)-2*I, " diff=", psi(0,z)-(log(z)-1/(2*z)-2*I))
# theta_{m}(z) := (psi(z) - [ln z - 1/(2z) - sum_{k<=m} B_{2k}/(2k z^{2k})]) / (-B_{2m+2}/((2m+2) z^{2m+2}))  should lie in (0,1)
def theta(m, z):
    A = log(z) - 1/(2*z) - sum(bernoulli(2*k)/(2*k*z**(2*k)) for k in range(1, m+1))
    return (psi(0, z) - A)/(-bernoulli(2*m+2)/((2*m+2)*z**(2*m+2)))
for m in [0,1,2,3,4,5,8]:
    ths = [theta(m, mpf(z)) for z in [1,2,3,5,10,100]]
    print(f"m={m}: theta(z) for z=1,2,3,5,10,100:", [float(t) for t in ths], " all in (0,1):", all(0 < t < 1 for t in ths))
# eps(d) = H_d - ln d - gamma = psi(d) + 1/d - ln d  -> eps(d) = 1/(2d) - sum B_2k/(2k d^2k) - theta B_{2m+2}/((2m+2) d^{2m+2})
for d in [1, 2, 7]:
    print("d=",d," eps(d)-(1/(2d)-sum_{k<=4}) =", harmonic(d)-log(d)-euler-(1/(2*mpf(d))-sum(bernoulli(2*k)/(2*k*mpf(d)**(2*k)) for k in range(1,5))), " theta_4(d)*(-B10/(10 d^10)) =", -theta(4,mpf(d))*bernoulli(10)/(10*mpf(d)**10))
# (b) non-squarefree: n with q=n/rad(n)>1: G(n)=(-1)^omega F(n) = sum_j c_j q^{-s_j} P_{s_j}(rad n) + E, |E|<= (1/132) q^{-10} prod(1+p^-10)
N = 3000
H = [Fraction(0)]
for i in range(1, N+1): H.append(H[-1] + Fraction(1, i))
def F_exact(n): return sum(int(mobius(n//d)) * H[d] for d in divisors(n))
def P(ps, s): return fprod([1 - mpf(p)**(-s) for p in ps])
worst = 0; cnt = 0; neg = []
for n in range(2, N+1):
    f = factorint(n); ps = sorted(f); k = len(ps)
    if k < 2: continue
    q = 1
    for p, e in f.items(): q *= p**(e-1)
    if q == 1: continue
    Gx = (-1)**k * F_exact(n); Gx = mpf(Gx.numerator)/Gx.denominator
    Gm = P(ps,1)/(2*q) - P(ps,2)/(12*q**2) + P(ps,4)/(120*q**4) - P(ps,6)/(252*q**6) + P(ps,8)/(240*q**8)
    Eb = fprod([1 + mpf(p)**(-10) for p in ps])/(132*mpf(q)**10)
    assert abs(Gx - Gm) <= Eb, (n, Gx, Gm, Eb)
    worst = max(worst, abs(Gx-Gm)/Eb); cnt += 1
    if Gx <= 0: neg.append(n)
print(f"non-squarefree n<={N}, omega>=2: {cnt} checked, formula holds, max |diff|/bound = {float(worst):.4f}; G<=0 for:", neg)
# (c) refine q=2 threshold: n = 2 * p_k#  (i.e. 4*3*5*...*p_k)
primes = list(primerange(2, 130000))
mp.dps = 30
def Pk(s, k): return fprod([1 - mpf(p)**(-s) for p in primes[:k]])
def G2(k):
    q = 2
    return Pk(1,k)/(2*q) - Pk(2,k)/(12*q**2) + Pk(4,k)/(120*q**4) - Pk(6,k)/(252*q**6) + Pk(8,k)/(240*q**8)
E2cap = mpf(1)/(132*1024)*(1+mpf('1e-3'))   # prod(1+p^-10) < zeta(10) < 1.001
lo, hi = 9250, 9500
while hi - lo > 1:
    mid = (lo+hi)//2
    if G2(mid) < 0: hi = mid
    else: lo = mid
print(f"q=2: G2({lo})={G2(lo)} > 0,  G2({hi})={G2(hi)} < 0; p_{hi}={primes[hi-1]}; error cap {E2cap}")
print("certified negative from k=", next(k for k in range(hi, hi+200) if G2(k) < -E2cap), " certified positive up to k=", max(k for k in range(lo-200, lo+1) if G2(k) > E2cap))
print("Q(k) at hi:", fprod([1+mpf(1)/p for p in primes[:hi]]))
