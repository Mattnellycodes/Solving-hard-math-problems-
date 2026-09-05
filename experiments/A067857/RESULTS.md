# A067857: Israel's sign conjecture is false in both directions — results summary

Date: 2026-09-05. Session: Claude Code (https://claude.ai/code/session_01SPbfjbe46J99BGSakrWEpq), workflow with
independent certificate / proof / audit sub-agents. All code and logs are under `experiments/A067857/`
(subdirectories listed in §6). Nothing here has been committed to git or posted anywhere.

## 0. One-paragraph summary

OEIS A067857 is defined by `sum_{d|n} a(d)/d! = H_n`, i.e. `a(n) = n! F(n)`, `F(n) = sum_{d|n} mu(n/d) H_d`.
Robert Israel's conjecture (OEIS comment, May 15 2015; `google-deepmind/formal-conjectures`,
`FormalConjectures/OEIS/67857.lean`, with `cardDistinctFactors` = omega = number of *distinct* prime factors) is

    a(n) < 0   <=>   omega(n) is odd and omega(n) >= 3.

**Both implications are false.** With `p_k#` the k-th primorial: `a(479#) < 0` although `omega(479#) = 92` is even
(`479 = p_92`, 199-digit n), and `a(487#) > 0` although `omega(487#) = 93` is odd (`487 = p_93`, 201-digit n).
More generally `sign a(p_k#) = (-1)^{k+1}` for **every** `k >= 92` — the exact opposite of the conjectured
`(-1)^k` — while the conjecture holds at every primorial with `k <= 91` and at every `n <= 411` (the OEIS b-file
range). These statements are (i) proved in `proof/PROOF.md` (Theorems A, B; the only numerical input is a finite
certificate for one inequality, `UB(92) < 0`, done in exact rational arithmetic with directed rounding), (ii)
certified numerically by three mutually independent implementations plus a standalone checker
(`certificate/`), and (iii) reproduced by two hostile audits with their own code and different methods
(`audit-certificate/`, `audit-proof/`), both of which returned verdict **CONFIRMED**.

## 1. Statement attacked (checked against the sources)

* OEIS entry fetched verbatim from the GitHub mirror `oeis/oeisdata` (`seq/A067/A067857.seq`): name
  "Sum_{k|n} a(k)/k! = Sum_{j=1 to n} 1/j"; data `1, 1, 5, 14, 154, 84, 8028, 25584, 361296, 528480, 80627040, 33471360`;
  comments "The first negative one is a(30) = -22690644647302814715858124800000" and
  "Conjecture: a(n) < 0 if and only if A001221(n) is an odd number >= 3. - Robert Israel, May 15 2015"
  (A001221 = omega).
* Lean formalisation (local clone of formal-conjectures, commit 8323e878):
  `a n := n! * sum_{d | n} mu(n/d) * harmonic d` and
  `theorem conjecture (n) (hn : 0 < n) : a n < 0 <-> Odd (cardDistinctFactors n) /\ 3 <= cardDistinctFactors n`
  (research open). This is exactly the statement refuted below. (A benchmark that formalised the statement with
  Omega, counted with multiplicity, and "refuted" it at n = 8 addressed a different, trivially false statement.)

## 2. What is PROVED (`proof/PROOF.md`, audited in `audit-proof/AUDIT.md`)

Notation: for squarefree n, `P_s(n) = prod_{p|n} (1 - p^{-s})`, `Q(n) = prod_{p|n} (1 + 1/p)`; for `omega(n) >= 2`,
`G(n) = (-1)^{omega(n)} F(n)`, so that the conjecture at squarefree n with omega >= 2 says exactly `G(n) > 0`
(for even omega it says `G(n) >= 0`; `a(n) = 0` never occurs in any certified range).

**Theorem A.** For every `k >= 92`, `sign a(p_k#) = (-1)^{k+1}`. In particular `a(479#) < 0` with omega even and
`a(487#) > 0` with omega odd >= 3; both directions of Israel's conjecture fail, infinitely often.

**Theorem B (corrected sign law for squarefree n, omega(n) = k >= 2, any m >= 1).**
`sign a(n) = (-1)^k sign G(n)` and

    G(n) = (1/2) P_1(n) - sum_{j=1}^m (B_{2j}/(2j)) P_{2j}(n) + c_m + E_m(n),
    c_m = 1/2 - gamma + sum_{j<=m} B_{2j}/(2j),
    |E_m(n)| <= (|B_{2m+2}|/(2m+2)) (prod_{p|n}(1 + p^{-2m-2}) - 1).

For m = 4: `G(n) = (1/2)P_1 - (1/12)P_2 + (1/120)P_4 - (1/252)P_6 + (1/240)P_8 + (2897/5040 - gamma) + E_4(n)`,
`|E_4(n)| < 1/132000`. Since `P_2 = P_1 Q`, this reads `G(n) = (1/2) P_1(n) (1 - Q(n)/6) + [0.0055 .. 0.0062] + E_4`,
which explains the phenomenon: the sign flips once `prod_{p|n}(1 + 1/p)` exceeds about 6.74 (for primorials:
`Q(91) = 6.7308 < 6.736 < Q(92) = 6.7448`).

Proof ingredients, all elementary except one: Moebius inversion; `sum_{d|n} mu(n/d) ln d = Lambda(n)` (so the
`ln d + gamma` parts cancel for omega >= 2); Euler products `sum_{d|n} mu(n/d) d^{-s} = (-1)^omega P_s(n)`;
`P_s(n) > 1/zeta(s)`; and the one analytic input, **Lemma 4** = DLMF 5.11(ii): for real z > 0 the Euler-Maclaurin
remainder of psi(z) is bounded by, and has the sign of, the first omitted term. PROOF.md re-proves Lemma 4 from
Binet's second formula for psi (so the DLMF citation is not load-bearing); the audit re-derived that proof by
hand and verified Binet's formula numerically to 1e-41 and theta in (0,1) for 150 (m, z) pairs at 150 digits.

Proof of Theorem A: with `UB(k) = (1/2)P_1(k) + (1/120)P_4(k) + (1/240)P_8(k) - 1/(2 pi^2) - 15/(4 pi^6) + 2897/5040 - gamma
+ 1/132000`, one has `G(p_k#) < UB(k)` (Lemma 5 applied to the negative terms), UB is strictly decreasing in k
(only P_1, P_4, P_8 depend on k), and **UB(92) <= -3.944e-5 < 0** by a finite certificate: exact rationals
`P_1(92), P_4(92), P_8(92)` over the 92 primes up to 479, `pi <= 3.1415926536` (upper bound), `gamma >= 0.5772156649`
(lower bound), eight terms each rounded up to 10 decimals, sum `-0.0000394407`. (Reference: `UB(91) = +5.49e-5`,
`UB(92) = -3.9440952182439e-5`, so the uniform bound singles out k = 92 exactly.) Hence `G(p_k#) < 0` for all
k >= 92. No Mertens-type theorem or any unproved asymptotic enters.

Further proved statements (Proposition 7 and Remark (d) of PROOF.md, thresholds certified by the audit):
* If `omega(n) >= 2`, `q = n / rad(n)` and `Q(rad n) <= 6q`, then `G(n) > 0`, i.e. the conjecture holds at n. Since
  `Q(51) = 5.98744 < 6 < Q(52) = 6.01249` and `Q(6480) = 11.99982 < 12 < Q(6481) = 12.0000032`: **every counterexample has
  omega(n) >= 52, and every non-squarefree counterexample has omega(n) >= 6481.** (Cases omega <= 1 are trivial:
  a(1) = 1, F(p^a) = H_{p^a} - H_{p^{a-1}} > 0.)
* Non-squarefree counterexamples exist too: for `n_k = 2 p_k# = 4*3*5*...*p_k`, `sign a(n_k) = (-1)^{k+1}` for all
  `k >= 9231` (`p_9231 = 95783`; n_9231 has 41,492 digits), and G(n_k) > 0 for 2 <= k <= 9230.

Audit verdict on the proof (`audit-proof/AUDIT.md`): CONFIRMED; every lemma re-derived by hand, the certificate for
UB(92) reproduced with MPFR directed rounding (gmpy2) and with hand-derived rational bounds for pi (Machin) and gamma
(Euler-Maclaurin at N = 1024 with the signed remainder), all remarks re-certified with interval arithmetic. One cosmetic
issue only: "iff G(n) > 0" should read "iff G(n) >= 0" when omega is even.

## 3. What is CERTIFIED numerically (`certificate/`, audited in `audit-certificate/README.md`)

### 3.1 Exact small n (n <= 411, the b-file range)
`a(n)` computed exactly in rational arithmetic both from the defining relation and by Moebius inversion; all integers,
agree, satisfy the defining relation, match the 20 OEIS data terms and a(30). **The conjecture holds for every
n <= 411** (87 negative indices, starting 30, 42, 60, 66, 70, ...; no zero terms). The Omega-with-multiplicity
variant fails first at n = 8. Reproduced by the audit with its own code (sympy + Fractions), and by the proof audit
for all n <= 4000 (conjecture holds; formula (B.2) error ratio < 1 at all 1882 squarefree n with omega >= 2).

### 3.2 Sign-definite enclosures of G(k) := (-1)^k F(p_k#) for 2 <= k <= 200
Method: closed form of Theorem B with m terms plus exact per-divisor corrections `eps(d) - A_m(d)` for all divisors
d <= D of p_k# (`eps(d) = H_d - ln d - gamma`), plus the rigorous tail `|E| <= |B_{2m+2}|/(2m+2) * sum_{d|n, d>D} d^{-2m-2}`
(an exact Euler-product identity). Transcendental inputs (ln p, gamma) are self-contained series enclosures with
explicit tail bounds; no library log/euler enters any certificate. Four implementations agree to all printed digits:

| implementation | arithmetic | (m, D) | widths |
|---|---|---|---|
| route (a) `route_a_iv.py` | mpmath.iv, 256-bit outward rounding | (4, 10^5) | ~3e-48 |
| route (b) `route_b_exact.py` | exact rationals only (H_d via lcm accumulator) | (6, 10^4) | ~3e-55 |
| checker `check_certificate.py` | pure-Python fixed-point (2^-400) integer intervals, re-derives everything from the JSON inputs | from JSON | passes, exit 0 |
| auditor's route C `audit-certificate/audit_routeC.py` | exact Fractions + own interval class on GNU MPFR (gmpy2, directed rounding); ln p by a telescoping atanh chain; gamma by Ein(256) - 8 ln 2 - E_1(256) (no Euler-Maclaurin at all) | (10, 3000) | ~1e-72, all 199 intervals strictly inside (a) and (b) |
| proof-auditor's Method D `audit-proof/audit_tight.py` | closed form + exact interval corrections for all d <= 60000, first-omitted-term tail | (4, 60000) | ~4e-39 |

Results: **G(k) > 0 for 2 <= k <= 91; G(k) < 0 for 92 <= k <= 200.** Critical values (route C, +-4.3e-73):

    G(91)  = +3.7316908332927197e-05   (a(467#) < 0, omega 91 odd: conjecture holds)
    G(92)  = -5.6769811107353809e-05   (a(479#) < 0, omega 92 even: 'only if' fails)
    G(93)  = -1.4912087756788696e-04   (a(487#) > 0, omega 93 odd:  'if' fails)
    G(94)  = -2.405327977e-04, G(100) = -7.592844578e-04, G(200) = -5.800641126e-03

Cross-checks that passed: exact rational brute force over all 2^k divisors for k <= 7 and 60-digit brute force for
k <= 14 lie inside every enclosure; F(6) = 7/60, F(30) = -199238355653/2329089562800 coincide with the exact part;
the checker detects all 8 tamperings injected by the auditor (flipped sign, shifted enclosure, tampered Euler product,
wrong prime, wrong Bernoulli number, wrong first-negative k, wrong status, shifted gamma), so it is not vacuous;
mpmath.iv basic operations and the decimal outward-rounding helpers were tested on 20000 random fractions.
Since Method D uses Euler-Maclaurin only for d > 60000, the sign of G(479#) does not depend on the constant in Lemma 4
in any delicate way.

### 3.3 Uniform bound, all k >= 92 (`certificate/uniform_bound.py`, `audit-certificate/audit_uniform_bound.py`)
Exact rationals only (zeta(2n) via Bernoulli numbers, pi via Machin resp. Euler's arctan formula, gamma enclosed).
For each m in {3,...,7}: UB_m strictly decreasing (checked to k = 300), UB_m(91) > 0 > UB_m(92)
(m = 4: +5.48e-5 / -3.95e-5), and `G(k) <= UB_m(k)` for every certified enclosure (min margin 7.8e-6). Together with
§3.2 this independently gives `sign a(p_k#) = (-1)^{k+1}` for all k >= 92.
(This is the same bound as Proposition 6 of the proof up to the choice of tail constant.)

## 4. What the audits found

* `audit-certificate` (hostile re-verification, own code, different method mix): **CONFIRMED.** All 199 enclosures
  reproduced strictly inside the certificate's; every derivation step checked; Lemma L1 re-proved from Binet's formula
  and checked at 260 digits; uniform bound reproduced; JSON bookkeeping cross-examined; checker mutation-tested.
  Negative result recorded by the auditor: a first run disagreed at k = 7 by 2.6e-18 — a gmpy2 unary-minus rounding
  bug in the auditor's own code, caught by the containment check; the certificate does not use gmpy2.
* `audit-proof`: **CONFIRMED.** Lemma-by-lemma table, Binet formula to 1e-41, theta in (0,1) at 150 digits, certificate
  table entries re-derived identically, all floating-point remarks of PROOF.md re-certified with directed rounding
  (margins >= 3e-8 against errors ~1e-40), (B.2)/(6.1) checked against exact rationals up to n = 4000.
* Residual reliance shared by everyone (not gaps in the argument): Binet's second formula for psi (DLMF 5.9;
  Whittaker-Watson 12.32), Euler's product and zeta(2n) formulas, `Ein(x) = gamma + ln x + E_1(x)` — textbook
  theorems, verified numerically to 40-110 digits, not re-proved from first principles. dlmf.nist.gov and oeis.org are
  egress-blocked in this sandbox: DLMF wording/equation numbers were confirmed via search snippets and independent
  derivation; the OEIS entry via the oeisdata GitHub mirror; the b-file itself was not compared (irrelevant, since
  part 3.1 verifies the defining relation exactly).

## 5. Open items and limits of what is claimed

1. **Minimality of 479# is not proved.** It is the first *primorial* counterexample (G(p_k#) > 0 for k <= 91 is
   certified), the conjecture holds for n <= 4000, and every counterexample has omega(n) >= 52 (Proposition 7); but a
   non-primorial squarefree n with 52 <= omega(n) <= 91 has not been excluded. The heuristic
   `G ~ (1/2)P_1(1 - Q/6) + 0.0055` and the fact that Q(n) is maximal at the primorial for fixed omega make a smaller
   counterexample look unlikely but this is not a proof. A natural next step is a proof that `G(n) >= G(p_omega(n)#)`
   or an explicit search over squarefree n with omega in [52, 91] and Q(n) > 6.7.
2. Full characterisation of the sign for general n (squarefree with all prime factors small, or non-squarefree with
   large omega) is only heuristic beyond Theorem B / (6.1); the corrected law "a(n) < 0 iff omega odd, omega >= 3,
   and G(n) > 0 (resp. omega even and G(n) < 0)" is a reformulation, not a closed criterion. A clean closed criterion
   would need control of Q(n) versus P_1(n), essentially `Q(rad n) <= 6 q` (sufficient, proved) versus something like
   `Q > 6.74 q` (necessary side, only heuristic).
3. The Lean statement in formal-conjectures is refuted mathematically, but no Lean proof of the negation was written;
   formalising Theorem A would require Lemma 4 (or Binet's formula) in Mathlib plus a certified evaluation of UB(92).
4. Nothing was posted to OEIS or committed; see the comment draft below.

## 6. Files

* `proof/PROOF.md` (Theorems A, B, Propositions 6, 7, certificate table), `proof/NOTES.md` (negative results:
  killed exact check for n up to 1e6, sharpness of the remainder bound, what the uniform bound cannot do),
  `proof/certificate.py`, `proof/verify1..5.py`.
* `certificate/README.md`, `route_a_iv.py`, `route_b_exact.py`, `check_certificate.py`, `cert_route_a.json`,
  `cert_route_b.json`, `part1_exact_small_n.py`, `uniform_bound.py`, `crosscheck_bruteforce.py`, `validate_lemmas.py`,
  `run_all.sh`, `logs/`.
* `audit-certificate/README.md`, `audit_routeC.py`, `results_routeC.json`, `audit_L1_binet.py`, `audit_uniform_bound.py`,
  `audit_small_n.py`, `cross_examine_json.py`, `toolchain_tests.py`, `mutation/`, `logs_*.txt`.
* `audit-proof/AUDIT.md`, `ivmpfr.py`, `audit_ub.py`, `audit_constants.py`, `audit_bruteforce.py`, `audit_lemma4.py`,
  `audit_lemma4b.py`, `audit_tight.py`, `*.out`.
* `main_agent_check/check.py` (the orchestrating agent's own 60-digit sanity check; not a certificate).

## 7. OEIS comment draft (<= 8 lines)

    The conjecture above is false in both directions: a(479#) < 0 although omega(479#) = 92 is even, and
    a(487#) > 0 although omega(487#) = 93 is odd (p_92 = 479, p_93 = 487); more precisely sign(a(p_k#)) = (-1)^(k+1)
    for all k >= 92, while the conjecture holds for all n <= 4000 and at every primorial p_k# with k <= 91.
    Reason: for squarefree n with omega(n) >= 2, a(n)/n! = (-1)^omega(n) * G(n) with
    G(n) = (1/2)*Product_{p|n}(1-1/p) - (1/12)*Product_{p|n}(1-1/p^2) + (1/120)*Product_{p|n}(1-1/p^4)
    - (1/252)*Product_{p|n}(1-1/p^6) + (1/240)*Product_{p|n}(1-1/p^8) + 2897/5040 - gamma + E, |E| < 1/132000
    (Euler-Maclaurin for H_d = psi(d)+1/d+gamma), so the sign flips once Product_{p|n}(1+1/p) exceeds about 6.74.
    Every counterexample has omega(n) >= 52. Computer-assisted (exact rationals + interval arithmetic), [name], Sep 05 2026.
