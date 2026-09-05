# Hostile audit of `experiments/A067857/proof/PROOF.md` (Theorems A, B; Propositions 6, 7; Remarks)

**Verdict: CONFIRMED.** Every lemma was re-derived by hand; the numeric certificate for Proposition 6(c) (UB(92) < 0)
was reproduced with a *different* method (MPFR directed rounding via gmpy2, plus hand-derived rational bounds for pi and
gamma), the signs for k = 91, 92, 93 were reproduced two ways, and no unjustified inequality was found. Minor remarks only
(section 4). All scripts in this directory are my own; outputs are the `*.out` files.

## 1. Sources examined
* `proof/PROOF.md`, `proof/NOTES.md`, `proof/verify1.py … verify5.py`, `proof/certificate.py`, `certificate.out`, `verify5.out`.
* Sibling work (read for cross-checking only): `certificate/README.md` + `certificate/logs/*`, `main_agent_check/check.py`.
* The formalisation, from the local clone of `google-deepmind/formal-conjectures` (commit 8323e878, 2026-09-04),
  file `FormalConjectures/OEIS/67857.lean`:
  `def a (n) : ℚ := if n = 0 then 0 else n! * ∑ d ∈ n.divisors, μ(n/d) * harmonic d` and
  `theorem conjecture (n) (hn : 0 < n) : a n < 0 ↔ Odd (cardDistinctFactors n) ∧ 3 ≤ cardDistinctFactors n` (research open).
  Docstring: "Conjecture: a(n) < 0 if and only if A001221(n) is an odd number ≥ 3" (A001221 = omega). This is exactly the
  statement attacked in PROOF.md (F(n) = a(n)/n! is the Lean definition itself, so Möbius inversion is not even needed).
* OEIS entry retrieved verbatim from the GitHub mirror `oeis/oeisdata` (`seq/A067/A067857.seq`): name "Sum_{k|n} a(k)/k! = Sum_{j=1 to n} 1/j"; data 1,1,5,14,154,84,8028,25584,361296,528480,80627040,33471360 (matches my exact computation, §3 C); comments "The first negative one is a(30) = -22690644647302814715858124800000" and "Conjecture: a(n) < 0 if and only if A001221(n) is an odd number >= 3. - Robert Israel, May 15 2015" (A001221 = omega). Statement and attribution in PROOF.md are correct.
* Scratch scripts in the scratchpad were read for orientation only.

## 2. Step-by-step check of the mathematics
| step | claim | check |
|---|---|---|
| §1 | F(n) = Σ μ(n/d) H_d, sign a = sign F | standard Möbius inversion; n! > 0. OK |
| Lemma 1(a) | Σ_{d\|n} μ(d) = 0, n > 1 | subsets of the ω(n) ≥ 1 primes, (1−1)^ω. OK |
| Lemma 1(b) | Σ μ(n/d) ln d = Λ(n) | ln n = Σ_{e\|n} Λ(e); swap sums; inner sum Σ_{g\|n/e} μ(g) = [e = n]. OK |
| Cor. 2 | F(n) = Σ μ(n/d) ε(d), ω(n) ≥ 2 | Λ(n) = 0, Σμ = 0 (n > 1). OK (also valid for non-squarefree n, used in Prop. 7) |
| Lemma 3 | Σ μ(n/d) f(d) = (−1)^k Σ (−1)^{ω(d)} f(d); Σ μ(n/d) d^{−s} = (−1)^k P_s(n) | μ(d_{S^c}) = (−1)^{k−\|S\|}; product expansion. OK (squarefree n only, as stated) |
| Lemma 4 / (3.3) | ψ(z) = ln z − 1/(2z) − 2∫ t/((t²+z²)(e^{2πt}−1)) dt | classical (differentiate Binet's 2nd formula: d/dz arctan(t/z) = −t/(t²+z²)); verified numerically to 1e−41 at z = 0.25, 0.5, 1, 2, 3.7, 10, 100 (`audit_lemma4.out`) |
| Lemma 4 | geometric expansion of 1/(t²+z²) with remainder (−1)^m t^{2m}/(z^{2m}(t²+z²)) | re-derived. OK |
| Lemma 4 | ∫ t^{2j+1}/(e^{2πt}−1) dt = (2j+1)! ζ(2j+2)/(2π)^{2j+2} = (−1)^j B_{2j+2}/(4(j+1)) > 0 | substitution x = 2πt + DLMF 25.5.1 + Euler's ζ(2n); verified numerically j = 0..5 to 40 digits |
| Lemma 4 | j-th term = −B_{2j+2}/((2j+2) z^{2j+2}); remainder = −θ B_{2m+2}/((2m+2) z^{2m+2}), θ ∈ (0,1) | algebra re-done; θ = ratio of integrals with 0 < z²/(t²+z²) < 1 pointwise. OK. Numerically θ_m(z) ∈ (0,1) for m ≤ 9, z ∈ {0.05,…,10^6} at 150 digits (`audit_lemma4b.out`); θ_4(1) = 0.3187, θ_4(2) = 0.6241, θ_4(71) = 0.99945 |
| (3.2) | ε(d) = ψ(d) + 1/d − ln d | H_d = ψ(d+1) + γ, ψ(d+1) = ψ(d) + 1/d. OK |
| Thm B | d = 1 term gives c_m; E_m over d > 1; \|E_m\| ≤ \|B_{2m+2}\|/(2m+2)·(Π(1+p^{−2m−2}) − 1) | R_m(1) = 1−γ−A_m(1) = c_m; Σ_{d\|n} d^{−s} = Π(1+p^{−s}). OK |
| Thm B, m = 4 | coefficients 1/12, −1/120, 1/252, −1/240; c_4 = 2897/5040 − γ; \|B_10\|/10 = 1/132 | 2520+420−42+20−21 = 2897 over 5040. OK |
| Thm B | ζ(10) − 1 ≤ 2^{−10} + 3^{−10} + 3^{−9}/9 = 0.00099914 < 10^{−3} | integral comparison Σ_{n≥4} n^{−10} ≤ ∫_3^∞. OK (ζ(10) − 1 = 0.00099458) |
| Lemma 5 | P_s(n) > 1/ζ(s), s > 1 | Euler product, missing factors (1−p^{−s})^{−1} > 1. OK |
| Prop. 6(a) | G < UB(k): (1/12)(6/π²) = 1/(2π²), (1/252)(945/π⁶) = 15/(4π⁶) | arithmetic OK; strict since P_2 > 6/π² |
| Prop. 6(b) | UB decreasing | P_1, P_4, P_8 multiplied by (1 − p^{−s}) ∈ (0,1). OK |
| Prop. 6(c) | UB(92) < 0 | certificate reproduced, see §3 |
| Thm A | sign a(p_k#) = (−1)^{k+1}, k ≥ 92 | G < UB(k) ≤ UB(92) < 0 ⇒ sign F = (−1)^k sign G = (−1)^{k+1}. OK; p_92 = 479, p_93 = 487 verified |
| Prop. 7 | squarefree: ½P_1(1−Q/6) ≥ 0, bracket ≥ 0.005459 | 90/(120π⁴) − 1/252 + 9450/(240π⁸) + c_4 − 1/132000 = 0.0054592 > 0. OK |
| Prop. 7 | non-squarefree: divisors with μ(n/d) ≠ 0 are d = qe, e \| rad n; (6.1); q^{−4}[0.007699 − 0.000992 − 0.000119] > 0 for q ≥ 2 | v_p(d) ∈ {a_p − 1, a_p}; no d = 1 term; P_6 ≤ 1, P_8 term dropped (≥ 0), Π(1+p^{−10}) < ζ(10) < 1.001. OK |
| Prop. 7 | Q(rad n) ≤ Q(ω(n)) | i-th prime of n ≥ p_i. OK. Thresholds certified in §3 |
| Remark (d) | R_4(2) kept exactly, residual < (1/132)2^{−10}(ζ(10)−1) < 7.4e−9; UB_2 decreasing | 7.36e−9. OK; certified in §3 |

## 3. Independent numerics (all mine; different tools from the proof)
**Method A – MPFR directed rounding (`ivmpfr.py`, `audit_ub.py`, 160 bits, every endpoint rounded outward; constants from
`const_pi`/`const_euler`/`log` with directed rounding).** The proof used exact `Fraction`s + a decimal table and `mpmath.iv`.
* p_91, p_92, p_93 = 467, 479, 487; exactly 92 primes ≤ 479.
* Exact rationals: P_1(92) = 0.0901580497409755666019668…, P_4(92) = 0.9239384033245916880001355…, P_8(92) = 0.9959392011255151468382712… (agree with (N2)).
* The eight table entries, re-derived as ceil(·10^10)/10^10 from the exact rationals and π ≤ 3.1415926536, γ ≥ 0.5772156649:
  identical to PROOF.md's row; sum = −0.0000394407 exactly as claimed. Exact-rational bound with those decimals: UB(92) ≤ −3.94409502e−5.
* UB(91) = +5.4866631374880875e−5, **UB(92) = −3.9440952182439463e−5**, UB(93) = −1.3200568511293e−4 (interval widths ≈ 1e−46);
  first k with UB(k) < 0 is 92; UB(k+1) < UB(k) numerically for 2 ≤ k < 200; G(p_k#) ≤ UB(k) for all k ≤ 200.
* Coarse enclosures with the per-n bound (B.1): G(p_k#) > 0 for all 2 ≤ k ≤ 91 (min lower bound 2.5071e−5 at k = 91);
  G(479#) ∈ [−6.9015e−5, −5.3960e−5]; G(487#) ∈ [−1.6137e−4, −1.4631e−4]; G < 0 for all 92 ≤ k ≤ 200. Matches Remark (a).
* Q(51) = 5.98744 < 6 < Q(52) = 6.01249; Q(6480) = 11.99982 < 12 < Q(6481) = 12.0000032 (p_6481 = 64853): certified, so
  Prop. 7's "ω ≥ 52" and "ω ≥ 6481" are rigorous.
* n_k = 2·p_k#: R_4(2) = −4.61753489087e−6; certified G(n_k) > 0 for 2 ≤ k ≤ 9230 (min lower bound 8.30e−8 at k = 9230),
  first certified G(n_k) < 0 at k = 9231 (p = 95783, upper bound −3.0e−8); UB_2(k) < 0 first at k = 9231. Matches Remark (d).

**Method B – hand-derived rational bounds for the constants (`audit_constants.py`).** Machin's formula with alternating-series
truncation gives π ∈ [3.14159265358979311599…, +9.6e−58]; Euler–Maclaurin at N = 1024, m = 4 with the *signed* remainder of
Lemma 4 and ln 2 = 2 atanh(1/3) gives γ ∈ [0.57721566490153286554…, +6e−33]. Hence π ≤ 3.1415926536 and γ ≥ 0.5772156649 (N3).

**Method C – exact rationals (`audit_bruteforce.py`).** a(n) from Σ a(d)/d! = H_n for n ≤ 40: integers, a(1..12) =
1, 1, 5, 14, 154, 84, 8028, 25584, 361296, 528480, 80627040, 33471360, a(30) = −22690644647302814715858124800000; Möbius form agrees.
(B.2) vs exact F(n) for all 1882 squarefree n ≤ 4000 with ω ≥ 2: max |error|/bound = 0.99924 (< 1); (6.1) for all 1528
non-squarefree n ≤ 4000 with ω ≥ 2: max ratio 0.99992 (< 1); the bound is essentially sharp, as the proof says. Conjecture holds
for every n ≤ 4000. Brute force over all 2^k divisors (50 digits) for k ≤ 12 agrees with the closed form within the bound.

**Method D – tight independent enclosure (`audit_tight.py`).** Closed form + exact interval corrections ε(d) − A_4(d) for all
divisors d ≤ 60000 of p_k# (interval H_d, MPFR log and γ) + first-omitted-term tail for d > 60000 (tail ≈ 1.4e−47). Widths 4e−39:
G(467#) = +3.7316908332927e−5, **G(479#) = −5.6769811107354e−5**, **G(487#) = −1.4912087756789e−4**, G(491#) = −2.4053e−4,
G(659#) = −2.2056e−3, G(1223#) = −5.8006e−3; exact rational values for k ≤ 6 lie inside. These agree to all printed digits with the
sibling `certificate/` enclosures (routes a/b). Note that here Euler–Maclaurin is only used for d > 60000, so the sign of G(479#)
would survive even a 10^40-fold weakening of Lemma 4's constant.

## 4. Issues found (none affects Theorems A/B)
1. **Cosmetic imprecision.** "The conjecture holds at a squarefree n with ω ≥ 2 iff G(n) > 0": for even ω the conjecture only
   predicts a(n) ≥ 0, so the correct statement is "iff G(n) ≥ 0" when ω is even (the difference is the case a(n) = 0, which never
   arises in the certified ranges).
2. **Citations.** dlmf.nist.gov is egress-blocked here too, so the DLMF wording/equation numbers could not be read verbatim.
   WebSearch snippets confirm: §5.11(ii) gives error bounds for the terminated expansions (5.11.1)/(5.11.2) (refs. Whittaker–Watson
   §12.33, Spira 1971, Olver 1997b pp. 293–295), and §5.9 contains the differentiated Binet formulas for ψ ((5.9.15)/(5.9.17) per the
   snippet), with the ψ-integral quoted verbatim in arXiv:1212.1432 as "frequently obtained by differentiating Binet's second
   formula". Since PROOF.md re-proves Lemma 4 from (3.3), the DLMF 5.11(ii) quotation is not load-bearing; (3.3) itself is a
   textbook formula (W&W §12.32 differentiated) and was verified numerically. The OEIS entry itself was retrieved verbatim from the
   `oeis/oeisdata` mirror (see §1); conjecture wording and attribution (Robert Israel, May 15 2015) are as stated in PROOF.md.
3. **Floating-point remarks now certified.** The author's caveat (2) — Remarks (a), (c) thresholds, (d) computed in plain
   30–40-digit floating point — is closed by Method A above (directed rounding, margins ≥ 3e−8 against errors ≈ 1e−40).
4. **Negative result of my own (recorded so the log is not misread).** `audit_lemma4.py` at 40 digits prints "FAIL" for
   θ_m(z) at z ≥ 500, m ≥ 4: pure cancellation (the remainder ≈ z^{−2m−2} < 10^{−40} is below working precision). At 150 digits
   (`audit_lemma4b.out`) all 150 (m, z) pairs give θ ∈ (0, 1), min 2.9e−4 (m = 9, z = 0.05), max 1 − 1e−13 (m = 0, z = 10^6).
5. Not proved (and not claimed): minimality of 479# among all n; only ω ≥ 52 / ω ≥ 6481 (non-squarefree) are established.

## 5. Files
`ivmpfr.py` (interval class), `audit_ub.py/.out`, `audit_constants.py/.out`, `audit_bruteforce.py/.out`, `audit_lemma4.py/.out`,
`audit_lemma4b.py/.out`, `audit_tight.py/.out`. Runtime: all under 25 s each.
