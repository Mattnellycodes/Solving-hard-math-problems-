#!/usr/bin/env python3
# Numerical audit of Lemma 4 / eq. (3.3): Binet-type integral for psi, the moment integrals, and theta_m(z) in (0,1).
from mpmath import mp, mpf, psi, log, quad, exp, pi, bernoulli, inf, zeta, gamma as Gam
mp.dps = 40
print("(3.3): psi(z) - [ln z - 1/(2z) - 2 int_0^inf t/((t^2+z^2)(e^{2 pi t}-1)) dt]")
for z in [mpf('0.25'), mpf('0.5'), mpf(1), mpf(2), mpf('3.7'), mpf(10), mpf(100)]:
    I = quad(lambda t: t/((t*t+z*z)*(exp(2*pi*t)-1)), [0, 1, 10, inf])
    print("  z=%6s  diff=%s" % (z, psi(0, z) - (log(z) - 1/(2*z) - 2*I)))
print("moment integrals int t^{2j+1}/(e^{2pi t}-1) dt vs (-1)^j B_{2j+2}/(4(j+1)):")
for j in range(0, 6):
    I = quad(lambda t: t**(2*j+1)/(exp(2*pi*t)-1), [0, 1, 10, inf])
    print("  j=%d  I=%s  formula=%s" % (j, I, (-1)**j*bernoulli(2*j+2)/(4*(j+1))))
def theta(m, z):
    A = log(z) - 1/(2*z) - sum(bernoulli(2*k)/(2*k*z**(2*k)) for k in range(1, m+1))
    return (psi(0, z) - A)/(-bernoulli(2*m+2)/((2*m+2)*z**(2*m+2)))
ok = True; worst = (0, None)
for m in range(0, 10):
    for z in [mpf('0.1'), mpf('0.3'), mpf('0.5'), mpf('0.9'), mpf(1), mpf('1.5'), mpf(2), mpf(3), mpf(5), mpf(10), mpf(71), mpf(500), mpf(10**4)]:
        th = theta(m, z)
        if not (0 < th < 1): ok = False; print("  FAIL m=%d z=%s theta=%s" % (m, z, th))
        if th > worst[0]: worst = (th, (m, z))
print("theta_m(z) in (0,1) for all tested (m,z):", ok, "; max theta:", worst)
print("theta_4(1) =", theta(4, mpf(1)), " theta_4(71) =", theta(4, mpf(71)), " theta_4(2) =", theta(4, mpf(2)))
print("Euler: zeta(2i) = (-1)^{i+1} B_{2i} (2pi)^{2i}/(2 (2i)!):", [zeta(2*i) - (-1)**(i+1)*bernoulli(2*i)*(2*pi)**(2*i)/(2*mp.factorial(2*i)) for i in range(1, 6)])
print("DLMF 25.5.1 check s=4: int x^3/(e^x-1) = Gamma(4) zeta(4):", quad(lambda x: x**3/(exp(x)-1), [0, inf]) - Gam(4)*zeta(4))
