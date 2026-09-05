# A067857 — rigorous sign certificates at primorials

OEIS A067857: `sum_{k|n} a(k)/k! = H_n`, equivalently `a(n) = n! F(n)`, `F(n) = sum_{d|n} mu(n/d) H_d`.
Conjecture (R. Israel, 2015; `FormalConjectures/OEIS/67857.lean`, with `cardDistinctFactors` = omega):

    a(n) < 0   <=>   omega(n) is odd and omega(n) >= 3.

## Results certified here

* **Part (1), exact.** For all n <= 411, a(n) computed exactly (rational arithmetic) from the defining relation and,
  independently, by Moebius inversion; both agree, satisfy the defining relation, are integers, match
  a(1..5) = 1, 1, 5, 14, 154 and a(30) = -22690644647302814715858124800000. The conjecture holds for every n <= 411
  (no a(n) = 0 occurs). The mis-formalised Omega-with-multiplicity variant fails already at n = 8 (102 failures below 411).
* **Part (2), rigorous enclosures of `G(k) := (-1)^k F(p_k#) = sum_{d | p_k#} (-1)^{omega(d)} H_d` for 2 <= k <= 200**
  by two independent routes and a third, standalone checker:
  * `G(k) > 0` for `2 <= k <= 91`  (conjecture holds at every primorial up to 467#),
  * `G(k) < 0` for `92 <= k <= 200`.
  Since `a(p_k#) = (p_k#)! (-1)^k G(k)`:
  * `n = 479# = p_92#` (omega = 92, even): **a(n) < 0** — violates the "only if" direction;
  * `n = 487# = p_93#` (omega = 93, odd): **a(n) > 0** — violates the "if" direction.
  Enclosure widths are ~3e-48 (route a) and ~3e-55 (route b), against |G(92)| = 5.68e-5, |G(91)| = 3.73e-5.
* **Bonus (uniform bound, exact rationals):** `G(k) <= UB(k) < 0` for **all** k >= 92 with UB strictly decreasing,
  so `sign a(p_k#) = (-1)^{k+1}` for every k >= 92: infinitely many counterexamples to both directions
  (`uniform_bound.py`, no Mertens-type theorem used; only pi and gamma enclosures and Euler products).

| k | p_k | omega parity | G(k) enclosure, route (a) [lo, hi] | width (a) | width (b) | sign G | sign a(p_k#) | conjecture predicts | status |
|---|-----|--------------|-------------------------------------|-----------|-----------|--------|--------------|---------------------|--------|
| 85 | 439 | odd | [+6.341513718e-04, +6.341513718e-04] | 3.0e-48 | 3.3e-55 | +1 | -1 | a<0 | conjecture holds at n=p_k# |
| 86 | 443 | even | [+5.310876427e-04, +5.310876427e-04] | 3.0e-48 | 3.3e-55 | +1 | +1 | a>=0 | conjecture holds at n=p_k# |
| 87 | 449 | odd | [+4.296278670e-04, +4.296278670e-04] | 3.0e-48 | 3.3e-55 | +1 | -1 | a<0 | conjecture holds at n=p_k# |
| 88 | 457 | even | [+3.301624306e-04, +3.301624306e-04] | 3.0e-48 | 3.3e-55 | +1 | +1 | a>=0 | conjecture holds at n=p_k# |
| 89 | 461 | odd | [+2.317742332e-04, +2.317742332e-04] | 3.1e-48 | 3.3e-55 | +1 | -1 | a<0 | conjecture holds at n=p_k# |
| 90 | 463 | even | [+1.340230285e-04, +1.340230285e-04] | 3.1e-48 | 3.3e-55 | +1 | +1 | a>=0 | conjecture holds at n=p_k# |
| 91 | 467 | odd | [+3.731690833e-05, +3.731690833e-05] | 3.1e-48 | 3.3e-55 | +1 | -1 | a<0 | conjecture holds at n=p_k# |
| 92 | 479 | even | [-5.676981111e-05, -5.676981111e-05] | 3.1e-48 | 3.3e-55 | -1 | -1 | a>=0 | VIOLATION of only if (a(n)<0 with omega even) |
| 93 | 487 | odd | [-1.491208776e-04, -1.491208776e-04] | 3.1e-48 | 3.2e-55 | -1 | +1 | a<0 | VIOLATION of if (a(n)>0 with omega odd >= 3) |
| 94 | 491 | even | [-2.405327977e-04, -2.405327977e-04] | 3.2e-48 | 3.2e-55 | -1 | -1 | a>=0 | VIOLATION of only if |
| 95 | 499 | odd | [-3.302989013e-04, -3.302989013e-04] | 3.2e-48 | 3.3e-55 | -1 | +1 | a<0 | VIOLATION of if |
| 96 | 503 | even | [-4.191738994e-04, -4.191738994e-04] | 3.2e-48 | 3.4e-55 | -1 | -1 | a>=0 | VIOLATION of only if |
| 97 | 509 | odd | [-5.068285884e-04, -5.068285884e-04] | 3.2e-48 | 3.4e-55 | -1 | +1 | a<0 | VIOLATION of if |
| 98 | 521 | even | [-5.923001429e-04, -5.923001429e-04] | 3.2e-48 | 3.5e-55 | -1 | -1 | a>=0 | VIOLATION of only if |
| 99 | 523 | odd | [-6.772817756e-04, -6.772817756e-04] | 3.3e-48 | 3.5e-55 | -1 | +1 | a<0 | VIOLATION of if |
| 100 | 541 | even | [-7.592844578e-04, -7.592844578e-04] | 3.3e-48 | 3.6e-55 | -1 | -1 | a>=0 | VIOLATION of only if |
| 125 | 691 | odd | [-2.521243966e-03, -2.521243966e-03] | 3.7e-48 | 4.5e-55 | -1 | +1 | a<0 | VIOLATION of if |
| 150 | 863 | even | [-3.858815976e-03, -3.858815976e-03] | 4.0e-48 | 5.7e-55 | -1 | -1 | a>=0 | VIOLATION of only if |
| 200 | 1223 | even | [-5.800641126e-03, -5.800641126e-03] | 4.5e-48 | 7.6e-55 | -1 | -1 | a>=0 | VIOLATION of only if |

(Full per-k data for k = 2..200 in `cert_route_a.json` / `cert_route_b.json`.)

## Mathematics (every inequality used)

Let n = p_k#, k >= 2, so n is squarefree, omega(n) = k, and mu(n/d) = (-1)^{k - omega(d)} for d | n. Hence
`F(n) = (-1)^k G(k)`, `G(k) = sum_{d|n} (-1)^{omega(d)} H_d`, and `sign a(n) = (-1)^k sign G(k)` (n! > 0).

**L1 (digamma remainder).** For real x > 0 and integer m >= 0,
`psi(x) = ln x - 1/(2x) - sum_{j=1}^m B_{2j}/(2j x^{2j}) + R_m(x)` where R_m(x) has the sign of `-B_{2m+2}` and
`|R_m(x)| < |B_{2m+2}| / ((2m+2) x^{2m+2})`. Reference: DLMF §5.11(ii) ("if z is real and positive, the remainder terms
in 5.11.2 are bounded in magnitude by the first neglected terms and have the same sign"); H. Alzer, *On some inequalities
for the gamma and psi functions*, Math. Comp. 66 (1997), Theorem 8 (the partial sums alternately over/under-estimate psi).
Since `H_d = psi(d) + 1/d + gamma`, with `eps(d) := H_d - ln d - gamma` and
`A_m(d) := 1/(2d) - sum_{j=1}^m B_{2j}/(2j d^{2j})` we get `eps(d) = A_m(d) + R_m(d)`.
(Checked numerically for 135 (x, m) pairs at 100 digits in `validate_lemmas.py`; this is the only non-elementary input.)

**L2.** For squarefree n with k >= 2 prime factors: `sum_{d|n} (-1)^{omega(d)} = (1-1)^k = 0`, and
`sum_{d|n} (-1)^{omega(d)} ln d = sum_{p|n} ln p * sum_{d: p | d | n} (-1)^{omega(d)} = sum_p ln p * (-(1-1)^{k-1}) = 0`.
Hence `G(k) = sum_{d|n} (-1)^{omega(d)} eps(d)`.

**L3.** By multiplicativity, `sum_{d|n} (-1)^{omega(d)} d^{-s} = prod_{p|n} (1 - p^{-s}) =: P_s` and
`sum_{d|n} d^{-s} = prod_{p|n} (1 + p^{-s})` (checked exactly for k <= 10, s <= 14).

**Closed form.** Splitting eps = A_m + R_m and the divisor set at a threshold D:

    G(k) = 1/2 P_1 - sum_{j=1}^m B_{2j}/(2j) P_{2j}                                (from A_m, L3)
           + sum_{d|n, d<=D} (-1)^{omega(d)} (eps(d) - A_m(d))                       (exact small-divisor part)
           + E,   |E| <= sum_{d|n, d>D} |R_m(d)| <= |B_{2m+2}|/(2m+2) * ( prod_{p|n}(1+p^{-(2m+2)}) - sum_{d|n, d<=D} d^{-(2m+2)} ).

The last bracket is exactly `sum_{d|n, d>D} d^{-(2m+2)}` (L3), so the bound is rigorous and tight.
For d <= D, `eps(d) - A_m(d)` is evaluated from H_d (exact rational or interval), ln d = sum_{p|d} ln p and gamma.

**Transcendental constants.** `ln p`: write p = 2^e y, y in [1,2), `ln y = 2 atanh(t)`, t = (y-1)/(y+1) in [0, 1/3],
`atanh(t) = sum_{j>=0} t^{2j+1}/(2j+1)`; the tail after J terms lies in `[0, t^{2J+1} / ((2J+1)(1-t^2))]`;
`ln 2 = 2 atanh(1/3)`. `gamma = H_D - ln D - A_m(D) - R_m(D)` with `|R_m(D)| < |B_{2m+2}|/((2m+2) D^{2m+2})` (L1 at x = D).
Both are self-contained rational/interval enclosures (widths 1e-97 and 1.5e-52 / 1.7e-57); mpmath's `iv.log` / `iv.euler`
are only used as consistency checks, never as inputs.

**Uniform bound (bonus).** With D = 1: `G(k) = sum_s c_s P_s(k) + (1 - gamma - A_m(1)) + E(k)`,
`|E(k)| <= C_m (prod_{i<=k}(1+p_i^{-(2m+2)}) - 1) <= C_m (zeta(2m+2)/zeta(4m+4) - 1)`. For s >= 2, `P_s(k)` decreases in k
and `P_s(k) >= 1/zeta(s)` (drop the factors of the convergent Euler product), so with S+ = {s: c_s > 0}, S- = {s: c_s < 0}:
`G(k) <= UB(k) := sum_{S+} c_s P_s(k) + sum_{S-} c_s/zeta(s) + (1 - gamma - A_m(1)) + C_m(zeta(2m+2)/zeta(4m+4) - 1)`,
strictly decreasing in k. `zeta(2n) = |B_{2n}|(2 pi)^{2n}/(2 (2n)!)`; pi is enclosed by Machin's formula (alternating series).
Result: UB(92) < 0 for each m in {3,...,7} (m = 4: UB(91) = +5.5e-5, UB(92) = -3.9e-5), hence G(k) < 0 for all k >= 92.

## Implementations (three independent kinds of arithmetic)

| file | arithmetic | (m, D) | what it does |
|---|---|---|---|
| `part1_exact_small_n.py` | exact fractions | – | a(n), F(n) for n <= 411 from the definition; conjecture check; a(30) |
| `route_a_iv.py` | mpmath.iv, 256-bit directed rounding | (4, 10^5) | Euler products for s = 1,2,4,6,8 as intervals; interval H_d (running sums); ln p by the series in interval arithmetic; gamma from H_D; tail bound over d > D; writes `cert_route_a.json` |
| `route_b_exact.py` | exact rationals only (no mpmath/sympy) | (6, 10^4) | H_d exact via an lcm accumulator; everything rational except the linear combination `sum_p c_p ln p + c_0 gamma` of series enclosures; writes `cert_route_b.json` |
| `check_certificate.py` | pure Python fixed-point integer intervals (2^-400) + fractions | from JSON | re-derives every enclosure from the JSON inputs (m, D, primes, Bernoulli numbers) with its own constants, verifies signs, consistency of all stored enclosures, prime list, divisor enumeration (brute-force count), Euler products |
| `crosscheck_bruteforce.py` | exact (k <= 7) / mpmath 60 digits (k <= 14) | – | brute force over all 2^k divisors; exact values must lie inside both enclosures |
| `validate_lemmas.py` | mpmath 100 digits + exact | – | numerical sanity checks of L1, L2, L3 and of the constant enclosures |
| `uniform_bound.py` | exact rationals | m = 3..7 | UB(k) < 0 for all k >= 92 |

Certificate JSON layout: `parameters`, `primes`, `bernoulli`, `coefficients`, `tail_constant`, `constants` (gamma, ln p enclosures
as outward-rounded decimal strings), `justification`, and `certificates[k]` with the Euler products, closed form, small-divisor
correction (route a) or the exact rational `Q(k)` with the integer coefficients `c_0`, `c_p` (route b), the tail bound `T`,
the final `G_enclosure`, `sign_G`, `sign_a`, `conjecture_predicts_a_negative`, `status`.

## Reproduction

    sh run_all.sh          # ~7 min total; route_b ~2 min, exact k=7 brute force ~3.5 min, everything else seconds
    python3 check_certificate.py cert_route_a.json cert_route_b.json   # standalone, ~8 s, exit code 0 = all verified

Logs of the runs used for this README are in `logs/`.

## Assumptions and caveats

* The only non-elementary ingredient is L1 (DLMF 5.11(ii) / Alzer 1997). Everything else is elementary and is either
  exact or done in outward-rounded interval arithmetic.
* Route (a) relies on mpmath.iv's basic operations (+, -, *, /, integer powers) being correctly rounded outward
  (they are implemented via `mpf_*` with `round_floor` / `round_ceiling`). No transcendental library function enters
  the certificate; `iv.log`/`iv.euler` appear only in assertions checking consistency. Route (b) and the checker use no
  floating point at all.
* Decimal strings in the JSON files are outward-rounded (floor for lower, ceiling for upper endpoints), so they are valid
  enclosures of the computed intervals; the checker parses them exactly (`fractions.Fraction`).
* Negative results / limits: the certificates cover primorials only (k <= 200 explicitly, all k >= 92 via the uniform bound).
  They say nothing about non-primorial n beyond n <= 411, where the conjecture holds. The OEIS b-file itself was not
  downloadable in this sandbox; part (1) is anchored on the published a(1..5) and a(30) and on the exact defining relation.
* The 60-digit mpmath brute force is a consistency check only (its own rounding is ~1e-59, visible as a 1e-59..1e-62
  "distance" from the essentially exact route-(b) enclosures for k <= 5); the rigorous brute-force containment test is the
  exact one for k <= 7.
