# Audit of `experiments/A067857/certificate` (hostile re-verification)

Auditor's verdict: **CONFIRMED** — every number, sign and inequality of the certificate was reproduced by
independent code with a different method mix; no unjustified inequality was found. The single non-elementary
input (lemma L1 = DLMF 5.11(ii)) is a standard theorem; it is re-derived below from Binet's second formula and
checked numerically at 260 digits. Residual reliance: Binet's second formula (DLMF §5.9, Whittaker–Watson §12.32)
and the classical constant evaluations used for the transcendental inputs (Euler's product for zeta,
A&S 5.1.11 for Ein/E_1) — textbook facts, checked numerically here to 60–110 digits.

## Files (all written by the auditor)

| file | what it does | result |
|---|---|---|
| `toolchain_tests.py` | mpmath.iv outward rounding, gmpy2 directed rounding, the certificate's `dec_floor/dec_ceil` on 20000 random (incl. negative) fractions | all pass |
| `audit_L1_binet.py`, `logs_L1_binet.txt`, `logs_L1_highprec.txt` | Binet's formula for psi (60 digits, 7 points), the moment integrals, L1 sign+magnitude for 247 (x,m) pairs (x from 0.1 to 10^5, m ≤ 12) at 110 digits, borderline cases redone at 260 digits, H_d = psi(d)+1/d+gamma, the Ein/E_1 identity | all pass (largest ratio bound 0.99999999913, sign always right; the 5 "FAIL"s at 110 digits were rounding noise, cleared at 260 digits) |
| `audit_small_n.py`, `results_small_n.json` | exact a(n), n ≤ 411, own code (sympy divisors/mobius + Fractions) | matches the 20 OEIS data terms (github.com/oeis/oeisdata), a(30), and every value/list in the certificate's `results_part1_small_n.json`; no violation for n ≤ 411; Omega-variant fails first at n = 8 |
| `audit_routeC.py`, `results_routeC.json`, `logs_routeC.txt` | **independent enclosure of G(k), k = 2..200**: exact Fractions for all algebraic parts; own interval class on GNU MPFR (gmpy2, RoundDown/RoundUp) for ln p and gamma; ln p by the chain ln p = Σ_{i=2}^{p} 2 atanh(1/(2i-1)); gamma by Ein(256) − 8 ln 2 − E_1(256) (no Euler–Maclaurin) and, as a cross-check, by L1 at 2^16; (m, D) = (10, 3000); divisors by factoring every integer ≤ D; brute-force recount for 8 values of k; brute force over all 2^k divisors for k ≤ 12 | widths ≈ 1e-72; **all 199 enclosures lie strictly inside both certificate enclosures**, signs agree; G>0 for k ≤ 91, G<0 for k ≥ 92 |
| `cross_examine_json.py`, `logs_cross_examine.txt` | the JSON files themselves: per-k divisor counts (spf-sieve recount), containment of exact P_s(k) in the stored Euler products, exact recomputation of their tail bounds T(k) with THEIR (m,D), bookkeeping of sign/parity/status, containment of my enclosures | 0 problems in both files; |G|/width ≥ 1.2e43 (a), 1.1e50 (b) |
| `mutation/` | the certificate's `check_certificate.py` on the original files (exit 0) and on 8 tampered copies (flipped sign, shifted enclosure, tampered Euler product, wrong prime, wrong B_14, wrong first-negative k, wrong status, shifted gamma) | all 8 tamperings detected (exit 1); the checker is not vacuous |
| `audit_uniform_bound.py`, `results_uniform_bound.json`, `logs_uniform_bound.txt` | own uniform bound UB_m(k) for m = 3..7, k ≤ 300: pi by Euler's arctan(1/2)+arctan(1/3) (not Machin), zeta(2n) by Bernoulli numbers, gamma via Ein; checks monotonicity, UB(92)<0<UB(91), G(k) ≤ UB(k) for my enclosures and for route a, and the Euler-product inequalities at k = 300 | reproduces the certificate's table to 4e-7 relative; m=4: UB(91)=+5.481830e-05, UB(92)=−3.948928e-05; strictly decreasing; G ≤ UB with margin ≥ 7.8e-6 |

Values at the critical indices (auditor's route C, outward-rounded; certificate values agree to all their printed digits):

    G(91) ∈ [+3.73169083329271974e-05 ± 4.3e-73]   (a(467#) < 0, omega = 91 odd: conjecture holds)
    G(92) ∈ [−5.67698111073538087e-05 ± 4.3e-73]   (a(479#) < 0, omega = 92 EVEN: 'only if' fails)
    G(93) ∈ [−1.491208775678869601e-04 ± 4.3e-73]  (a(487#) > 0, omega = 93 odd ≥ 3: 'if' fails)

## The mathematics, re-derived

Notation: n = p_k#, k ≥ 2; F(n) = Σ_{d|n} μ(n/d) H_d = a(n)/n!; μ(n/d) = (−1)^{k−ω(d)}, so F(n) = (−1)^k G(k),
G(k) = Σ_{d|n} (−1)^{ω(d)} H_d, and sign a(n) = (−1)^k sign G(k). (Checked against exact a(6), a(30), a(210).)

*L2.* Σ_{d|n}(−1)^{ω(d)} = (1−1)^k = 0 and Σ_{d|n}(−1)^{ω(d)} ln d = Σ_{p|n} ln p · Σ_{d': d'|n/p} (−1)^{1+ω(d')} = −(1−1)^{k−1} Σ ln p = 0
for k ≥ 2. Hence G(k) = Σ (−1)^{ω(d)} ε(d), ε(d) = H_d − ln d − γ.

*L3.* Σ_{d|n}(−1)^{ω(d)} d^{−s} = Π_{p|n}(1 − p^{−s}), Σ_{d|n} d^{−s} = Π_{p|n}(1 + p^{−s}) (multiplicativity; verified exactly for k ≤ 12).

*L1 (own proof).* Binet's second formula for ψ, x > 0 (DLMF §5.9; obtained by differentiating Binet's second formula for ln Γ
under the integral sign, Whittaker–Watson §12.32; verified here numerically to 1e-60 at 7 points):
ψ(x) = ln x − 1/(2x) − 2∫_0^∞ t dt / ((t²+x²)(e^{2πt} − 1)).
Insert 1/(t²+x²) = Σ_{j=0}^{m−1} (−1)^j t^{2j}/x^{2j+2} + (−1)^m t^{2m}/(x^{2m}(t²+x²)) and use
∫_0^∞ t^{2j+1}/(e^{2πt}−1) dt = (2j+1)! ζ(2j+2)/(2π)^{2j+2} = |B_{2j+2}|/(4(j+1)) = (−1)^j B_{2j+2}/(4(j+1))
(expand 1/(e^{2πt}−1) = Σ_{r≥1} e^{−2πrt}, integrate termwise — monotone convergence — and use Euler's ζ(2n) formula;
verified numerically to 1e-59 for j ≤ 7). The j-th term contributes −B_{2j+2}/((2j+2) x^{2j+2}), i.e. the Euler–Maclaurin
sum, and the remainder is R_m(x) = −2(−1)^m x^{−2m} ∫_0^∞ t^{2m+1} dt/((t²+x²)(e^{2πt}−1)). Because 0 < 1/(t²+x²) < 1/x² on t > 0
and the rest of the integrand is positive, R_m(x) = −θ B_{2m+2}/((2m+2) x^{2m+2}) with θ ∈ (0,1): R_m has the sign of −B_{2m+2}
and |R_m(x)| < |B_{2m+2}|/((2m+2) x^{2m+2}). This is exactly DLMF 5.11(ii) ("when z is real and positive, the remainder terms
are bounded in magnitude by the first neglected terms and have the same sign" — wording confirmed by a web-search snippet of
dlmf.nist.gov/5.11; the site itself is blocked from this sandbox). With H_d = ψ(d) + 1/d + γ: ε(d) = A_m(d) + R_m(d),
A_m(d) = 1/(2d) − Σ_{j≤m} B_{2j}/(2j d^{2j}). Numerically (260 digits) the ratio |R_m|/bound is < 1 with the right sign for all
tested (x, m), tending to 1 from below as x → ∞, as the proof predicts.

*Closed form.* G(k) = Σ_s c_s P_s(k) + Σ_{d≤D, d|n} (−1)^{ω(d)} (ε(d) − A_m(d)) + E, c_1 = 1/2, c_{2j} = −B_{2j}/(2j),
|E| ≤ C_m Σ_{d|n, d>D} d^{−(2m+2)} = C_m (Π_{p|n}(1+p^{−(2m+2)}) − Σ_{d≤D,d|n} d^{−(2m+2)}), C_m = |B_{2m+2}|/(2m+2).
For d ≤ D, ε(d) − A_m(d) = (H_d − A_m(d)) − Σ_{p|d} ln p − γ, so only integer multiples of ln p and γ are transcendental.
Verified: with (m, D) = (10, 3000) the closed form reproduces brute-force G(k) for k ≤ 12 to 1e-50 and, in mpmath, G(7) to 2e-59
independently of D ∈ {3000, 3100, 10^4, 6·10^5}.

*Transcendental inputs (auditor's, different from the certificate's).* ln p = Σ_{i=2}^{p} 2 atanh(1/(2i−1)), atanh(1/n) =
Σ_j 1/((2j+1)n^{2j+1}) with tail ≤ first omitted term · n²/(n²−1); cross-checked against MPFR's correctly-rounded log and against
the certificate's 200 enclosures (widths 4e-123 vs 1.4e-97). γ = Ein(256) − 8 ln 2 − E_1(256), Ein(x) = Σ_{k≥1} (−1)^{k−1}x^k/(k·k!)
(exact alternating series, remainder ≤ first omitted term once k > x), 0 < E_1(x) = ∫_x^∞ e^{−t}/t dt < e^{−x}/x < 2.7^{−x}/x
(e > 2.7); identity Ein(x) = γ + ln x + E_1(x) is A&S 5.1.11 / DLMF 6.2.3–6.2.4 (derivative of the difference is 0, constant fixed by
Γ'(1) = −γ; verified numerically to 1e-111). This γ (width 1.5e-113) agrees with the L1-based γ at N = 2^16 (width 6e-104), with
the certificate's γ (from H_{10^5}, H_{10^4}) and with the published digits — an independent confirmation of L1 at x = 10^4, 10^5, 2^16.

*Uniform bound.* With D = 1: G(k) = Σ c_s P_s(k) + (1 − γ − A_m(1)) + E(k), |E(k)| ≤ C_m(Π_{i≤k}(1+p_i^{−(2m+2)}) − 1)
≤ C_m(ζ(2m+2)/ζ(4m+4) − 1). For s > 1, P_s(k) decreases in k and P_s(k) ≥ Π_all(1−p^{−s}) = 1/ζ(s); Π_{i≤k}(1+p_i^{−s}) ≤ ζ(s)/ζ(2s).
So G(k) ≤ UB_m(k) = Σ_{c_s>0} c_s P_s(k) + Σ_{c_s<0} c_s/ζ(s) + (1 − γ − A_m(1)) + C_m(ζ(2m+2)/ζ(4m+4) − 1), strictly decreasing
(c_1 = 1/2 > 0, P_1 strictly decreasing). Directions of rounding in the certificate's code checked line by line: c_s<0 uses the
upper bound of ζ(s); γ enters with its lower bound; the tail uses upper ζ(2m+2) / lower ζ(4m+4). Reproduced: UB_m(92) < 0 for
m = 3..7, monotone to k = 300, hence G(k) < 0 and sign a(p_k#) = (−1)^{k+1} for every k ≥ 92 (no Mertens-type input).

## Negative results and limitations (recorded)

* My first route-C run disagreed with the certificate at k = 7 by 2.6e-18. The bug was **mine**: gmpy2 evaluates Python unary
  minus on an mpfr in the global 53-bit context, so `-x` silently rounded a 420-bit interval endpoint. It surfaced only when a
  coefficient c_p < 0 first appears (k = 7) and was caught by the containment check against the certificate (the certificate does
  not use gmpy2; the bug does not affect it). Fixed by multiplying by the signed exact coefficient under the directed contexts.
* Two further failures were artifacts of my sanity assertions (truncated 70-digit γ literal below a 1e-113-wide enclosure; 60-digit
  output strings wider than the certificate's 80-digit ones at k ≤ 5); fixed, no effect on conclusions.
* dlmf.nist.gov, oeis.org and ams.org are blocked by the egress proxy: the DLMF 5.11(ii) wording is confirmed only through a search
  snippet (and by the proof above); Alzer 1997 "Theorem 8" (a secondary citation in the certificate) could not be checked and is
  not load-bearing; the OEIS b-file (n ≤ 411) was not available — only the 20 data terms and a(30) were compared, but part (1) is
  anchored on the defining relation itself, so the b-file is irrelevant to correctness.
* Scope: the certificate (and this audit) covers all n ≤ 411 and the primorials p_k#, 2 ≤ k ≤ 200 explicitly plus all k ≥ 92 by the
  uniform bound. Nothing is claimed or audited about other n > 411.
