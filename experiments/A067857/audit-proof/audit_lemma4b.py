#!/usr/bin/env python3
# Re-check theta_m(z) in (0,1) at 150 digits (the 40-digit run showed cancellation artefacts for z >= 500, where the
# remainder ~ z^{-2m-2} is below working precision).
from mpmath import mp, mpf, psi, log, bernoulli
mp.dps = 150
def theta(m, z):
    A = log(z) - 1/(2*z) - sum(bernoulli(2*k)/(2*k*z**(2*k)) for k in range(1, m+1))
    return (psi(0, z) - A)/(-bernoulli(2*m+2)/((2*m+2)*z**(2*m+2)))
ok = True; mx = (0, None); mn = (1, None)
zs = [mpf('0.05'), mpf('0.1'), mpf('0.3'), mpf('0.5'), mpf('0.9'), mpf(1), mpf('1.5'), mpf(2), mpf(3), mpf(5), mpf(10), mpf(71), mpf(500), mpf(10**4), mpf(10**6)]
for m in range(0, 10):
    for z in zs:
        th = theta(m, z)
        if not (0 < th < 1): ok = False; print("  FAIL m=%d z=%s theta=%s" % (m, z, mp.nstr(th, 20)))
        if th > mx[0]: mx = (th, (m, z))
        if th < mn[0]: mn = (th, (m, z))
print("theta_m(z) in (0,1) for all m<=9 and all tested z (150 digits):", ok)
print("max theta:", mp.nstr(mx[0], 15), mx[1], " min theta:", mp.nstr(mn[0], 15), mn[1])
print("theta_4(d) for d=1..8:", [mp.nstr(theta(4, mpf(d)), 12) for d in range(1, 9)])
