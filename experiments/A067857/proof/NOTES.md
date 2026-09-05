# Working notes for experiments/A067857/proof (negative results and verification log)

Scripts (all mine, written for this task): `verify1.py` (killed: exact-rational check of F(n) for n up to ~10^6 via
Fraction harmonic numbers is far too slow; it did confirm a(30) = OEIS value and no violations for n<200 before being killed),
`verify1b.py` (exact-rational check of (B.2) for all 2862 squarefree n<=6000 with omega>=2, and 40-digit brute force over all
2^k divisors for primorials k<=12), `verify2.py` (UB(k), certified enclosures of G(p_k#), k<=91 positivity, monotonicity, heuristic
constants), `verify3.py` (Binet representation to 1e-41, theta in (0,1), non-squarefree formula (6.1) vs exact rationals for 1140 n<=3000,
q=2 threshold), `verify4.py` (interval arithmetic for UB(92), coarse 10-decimal certificate, k6=52, k12=6481, Prop. 7 constants,
n_k=2 p_k# flip at k=9231; its last O(k^2) loop timed out and was redone in `verify5.py` incrementally), `certificate.py`
(exact-rational certificate). Outputs: `certificate.out`, `verify5.out`.

Negative results / limitations
- The uniform bound UB(k) cannot certify anything for k<=91 (UB(91)=+5.49e-5>0); positivity for k<=91 needs the per-n error bound
  and P_s(k) for s in {1,2,4,6,8,10} (done: min margin 2.5e-5 at k=91).
- The first-omitted-term bound |R_4(d)| <= 1/(132 d^10) is essentially sharp (ratio 0.9995 attained at n=71*73, theta_4(71)~1), so the
  1/132000 cap cannot be improved by better bounding of the same remainder; only by larger m or by treating small d exactly.
- For n_k = 2 p_k#, with the 1/(132*1024)*1.001 cap alone the flip could only be located to k in [9209, 9325]; keeping R_4(2) exactly
  (error cap 7.4e-9) pins it to k=9231 and the uniform bound UB_2 also turns negative exactly at 9231.
- dlmf.nist.gov is blocked by the egress proxy; equation numbers 5.11.2/5.11(ii), 5.5.2, 25.2.11, 25.5.1, 25.6.2 are quoted from
  memory and search snippets (5.11(ii) wording and 25.5.1 confirmed by search snippets). The psi Binet-type integral is cited by
  section (DLMF §5.9) and Whittaker–Watson §12.32, because its DLMF equation number (5.9.15 or 5.9.16) could not be confirmed.
- mpmath.iv does not work with fprod (raised ValueError); products were formed with an explicit loop.
- Minimality of 479# among ALL counterexamples is not proved (only: every counterexample has omega>=52, non-squarefree ones omega>=6481,
  and 479# is the first primorial counterexample).
