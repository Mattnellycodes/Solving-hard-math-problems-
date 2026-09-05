# A067857 at primorials: the sign of $a(p_k\#)$ and the failure of Israel's conjecture

**Conventions.** $p_1=2<p_2=3<\cdots$ are the primes and $p_k\#=p_1p_2\cdots p_k$. $H_n=\sum_{j\le n}1/j$, $\gamma$ is Euler's constant, $\psi=\Gamma'/\Gamma$, $B_j$ are the Bernoulli numbers ($B_2=\tfrac16,\ B_4=-\tfrac1{30},\ B_6=\tfrac1{42},\ B_8=-\tfrac1{30},\ B_{10}=\tfrac5{66}$), $\mu$ is the Möbius function, $\Lambda$ the von Mangoldt function, $\omega(n)$ the number of *distinct* prime factors, $\operatorname{rad}(n)=\prod_{p\mid n}p$.
OEIS A067857: $a(n)$ is defined by $\sum_{d\mid n}a(d)/d!=H_n$. **Israel's conjecture** (OEIS comment 2015; `formal-conjectures`, `OEIS/67857.lean`, with `cardDistinctFactors` $=\omega$): $a(n)<0 \iff \omega(n)$ is odd and $\omega(n)\ge 3$. (The variant with $\Omega$ in place of $\omega$ is a different statement, false already at $n=8$; it plays no role here.)

For squarefree $n$ and real $s\ge1$ put $P_s(n)=\prod_{p\mid n}(1-p^{-s})$ and $Q(n)=\prod_{p\mid n}(1+1/p)$, so that $P_2(n)=P_1(n)Q(n)$; abbreviate $P_s(k)=P_s(p_k\#)$, $Q(k)=Q(p_k\#)$. Put $F(n)=a(n)/n!$ and, for $\omega(n)\ge2$, $G(n)=(-1)^{\omega(n)}F(n)$.

**Theorem A.** For every $k\ge 92$, $\operatorname{sign}a(p_k\#)=(-1)^{k+1}$. In particular $a(479\#)<0$ although $\omega(479\#)=92$ is even, and $a(487\#)>0$ although $\omega(487\#)=93$ is odd and $\ge3$ ($p_{92}=479$, $p_{93}=487$). Both implications of the conjecture are false.

**Theorem B (corrected sign law).** Let $n$ be squarefree with $\omega(n)=k\ge2$ and let $m\ge1$. Then $\operatorname{sign}a(n)=(-1)^k\operatorname{sign}G(n)$ and
$$G(n)=\tfrac12P_1(n)-\sum_{j=1}^m\frac{B_{2j}}{2j}P_{2j}(n)+c_m+E_m(n),\qquad c_m=\tfrac12-\gamma+\sum_{j=1}^m\frac{B_{2j}}{2j},\qquad |E_m(n)|\le\frac{|B_{2m+2}|}{2m+2}\Big(\prod_{p\mid n}(1+p^{-2m-2})-1\Big),\tag{B.1}$$
where $E_m(n)=\sum_{d\mid n,\,d>1}(-1)^{\omega(d)}R_m(d)$ with $R_m(d)$ the remainder of the Euler–Maclaurin expansion (3.2) of $H_d-\ln d-\gamma$.
For $m=4$:
$$G(n)=\tfrac12P_1-\tfrac1{12}P_2+\tfrac1{120}P_4-\tfrac1{252}P_6+\tfrac1{240}P_8+c_4+E_4(n),\quad c_4=\tfrac{2897}{5040}-\gamma=-0.0024140776\ldots,\quad |E_4(n)|<\tfrac1{132000}<7.6\cdot10^{-6}.\tag{B.2}$$
Since the conjecture predicts $\operatorname{sign}a(n)=(-1)^{k}$ for such $n$, it holds at a squarefree $n$ with $\omega(n)\ge2$ **iff** $G(n)>0$.

## 1. Möbius inversion and elimination of $\ln d+\gamma$

Möbius inversion of $\sum_{d\mid n}a(d)/d!=H_n$ gives $F(n)=a(n)/n!=\sum_{d\mid n}\mu(n/d)H_d$; in particular $\operatorname{sign}a(n)=\operatorname{sign}F(n)$.

**Lemma 1.** (a) $\sum_{d\mid n}\mu(d)=0$ for $n>1$. (b) $\sum_{d\mid n}\mu(n/d)\ln d=\Lambda(n)$ for all $n\ge1$; in particular the sum vanishes when $\omega(n)\ge2$.

*Proof.* (a) The divisors $d$ of $n$ with $\mu(d)\ne0$ are the products $d_S=\prod_{p\in S}p$ over subsets $S$ of the $\omega(n)$ primes dividing $n$, and $\mu(d_S)=(-1)^{|S|}$; so the sum is $(1-1)^{\omega(n)}=0$.
(b) If $n=\prod_p p^{a_p}$ then $\ln n=\sum_p a_p\ln p=\sum_{p}\sum_{1\le j\le a_p}\Lambda(p^j)=\sum_{e\mid n}\Lambda(e)$. Hence
$\sum_{d\mid n}\mu(n/d)\ln d=\sum_{d\mid n}\mu(n/d)\sum_{e\mid d}\Lambda(e)=\sum_{e\mid n}\Lambda(e)\sum_{e\mid d\mid n}\mu(n/d)$. For fixed $e\mid n$, as $d$ runs through the multiples of $e$ dividing $n$, $n/d$ runs through the divisors of $n/e$; so the inner sum is $\sum_{g\mid n/e}\mu(g)$, which by (a) is $1$ if $e=n$ and $0$ otherwise. The total is $\Lambda(n)$, which is $0$ unless $n$ is a prime power. $\square$

**Corollary 2.** If $\omega(n)\ge2$ then $F(n)=\sum_{d\mid n}\mu(n/d)\,\varepsilon(d)$, where $\varepsilon(d):=H_d-\ln d-\gamma$.
*Proof.* $\sum_{d\mid n}\mu(n/d)(\ln d+\gamma)=\Lambda(n)+\gamma\sum_{d\mid n}\mu(d)=0+0$. $\square$

## 2. Euler products

**Lemma 3.** Let $n$ be squarefree with $\omega(n)=k$ and let $f$ be any function on the divisors of $n$. Then $\sum_{d\mid n}\mu(n/d)f(d)=(-1)^k\sum_{d\mid n}(-1)^{\omega(d)}f(d)$, and for $f(d)=d^{-s}$ ($s>0$): $\ \sum_{d\mid n}\mu(n/d)\,d^{-s}=(-1)^kP_s(n)$.

*Proof.* Divisors of $n$ correspond to subsets $S$ of its $k$ primes, $d=d_S$; then $n/d_S=d_{S^c}$ and $\mu(n/d_S)=(-1)^{k-|S|}=(-1)^k(-1)^{\omega(d_S)}$. For $f(d)=d^{-s}$: $\sum_S(-1)^{|S|}\prod_{p\in S}p^{-s}=\prod_{p\mid n}(1-p^{-s})$. $\square$

## 3. The expansion of $\varepsilon(d)$ with a signed remainder

**Lemma 4.** For real $z>0$ and integer $m\ge0$ there is $\theta=\theta_m(z)\in(0,1)$ with
$$\psi(z)=\ln z-\frac1{2z}-\sum_{j=1}^m\frac{B_{2j}}{2j\,z^{2j}}-\theta\,\frac{B_{2m+2}}{(2m+2)\,z^{2m+2}}.\tag{3.1}$$
Consequently, for every integer $d\ge1$,
$$\varepsilon(d)=\frac1{2d}-\sum_{j=1}^m\frac{B_{2j}}{2j\,d^{2j}}+R_m(d),\qquad R_m(d)=-\theta_m(d)\frac{B_{2m+2}}{(2m+2)d^{2m+2}},\qquad |R_m(d)|\le\frac{|B_{2m+2}|}{(2m+2)\,d^{2m+2}}.\tag{3.2}$$

(3.1) with $\theta\in[0,1]$ is exactly the statement of **DLMF §5.11(ii)** for the expansion DLMF 5.11.2, $\psi(z)\sim\ln z-\frac1{2z}-\sum_{k\ge1}B_{2k}/(2kz^{2k})$: *"If the sums in the expansions (5.11.1) and (5.11.2) are terminated at $k=n-1$ ($k\ge0$) and $z$ is real and positive, then the remainder terms are bounded in magnitude by the first neglected terms and have the same sign."* (3.2) follows from (3.1) because $H_d=\psi(d+1)+\gamma=\psi(d)+1/d+\gamma$ (DLMF 5.4.14 and 5.5.2), so $\varepsilon(d)=\psi(d)-\ln d+1/d$.

*Proof of (3.1) from Binet's second formula* (for completeness). For $z>0$, differentiating Binet's second formula for $\ln\Gamma$ under the integral sign (DLMF §5.9; Whittaker–Watson §12.32) gives
$$\psi(z)=\ln z-\frac1{2z}-2\int_0^\infty\frac{t\,dt}{(t^2+z^2)(e^{2\pi t}-1)}.\tag{3.3}$$
For $t>0$, $\frac1{t^2+z^2}=\sum_{j=0}^{m-1}\frac{(-1)^jt^{2j}}{z^{2j+2}}+\frac{(-1)^mt^{2m}}{z^{2m}(t^2+z^2)}$. By DLMF 25.5.1, $\int_0^\infty x^{s-1}(e^x-1)^{-1}dx=\Gamma(s)\zeta(s)$, and Euler's formula $\zeta(2i)=(-1)^{i+1}B_{2i}(2\pi)^{2i}/(2(2i)!)$ (DLMF 25.6.2),
$$\int_0^\infty\frac{t^{2j+1}dt}{e^{2\pi t}-1}=\frac{(2j+1)!\,\zeta(2j+2)}{(2\pi)^{2j+2}}=\frac{(-1)^jB_{2j+2}}{4(j+1)}>0 .$$
Inserting the expansion into (3.3), the $j$-th term contributes $-2(-1)^jz^{-2j-2}\cdot(-1)^jB_{2j+2}/(4(j+1))=-B_{2j+2}/((2j+2)z^{2j+2})$, i.e. the sum in (3.1), and the remainder is $-2(-1)^mz^{-2m}\int_0^\infty t^{2m+1}\big((t^2+z^2)(e^{2\pi t}-1)\big)^{-1}dt$. Since $0<(t^2+z^2)^{-1}<z^{-2}$ for $t>0$ and the remaining factor of the integrand is positive, this equals $\theta\cdot\big(-2(-1)^mz^{-2m-2}(-1)^mB_{2m+2}/(4(m+1))\big)=-\theta B_{2m+2}/((2m+2)z^{2m+2})$ with $\theta\in(0,1)$. $\square$

## 4. Proof of Theorem B

Let $n$ be squarefree with $\omega(n)=k\ge2$. By Corollary 2 and (3.2),
$F(n)=\sum_{d\mid n}\mu(n/d)\big[\tfrac1{2d}-\sum_{j\le m}\tfrac{B_{2j}}{2j d^{2j}}\big]+\sum_{d\mid n}\mu(n/d)R_m(d)$.
By Lemma 3 the first sum is $(-1)^k\big[\tfrac12P_1(n)-\sum_{j\le m}\tfrac{B_{2j}}{2j}P_{2j}(n)\big]$ and the second is $(-1)^k\sum_{d\mid n}(-1)^{\omega(d)}R_m(d)$. Its term $d=1$ is $(-1)^kR_m(1)=(-1)^k\big(\varepsilon(1)-\tfrac12+\sum_{j\le m}\tfrac{B_{2j}}{2j}\big)=(-1)^kc_m$ because $\varepsilon(1)=1-\gamma$. Put $E_m(n):=\sum_{d\mid n,\,d>1}(-1)^{\omega(d)}R_m(d)$; by (3.2), $|E_m(n)|\le\frac{|B_{2m+2}|}{2m+2}\sum_{d\mid n,d>1}d^{-2m-2}=\frac{|B_{2m+2}|}{2m+2}\big(\prod_{p\mid n}(1+p^{-2m-2})-1\big)$. Multiplying by $(-1)^k$ gives (B.1), and $\operatorname{sign}a(n)=\operatorname{sign}F(n)=(-1)^k\operatorname{sign}G(n)$.
For $m=4$: $B_2/2=\tfrac1{12}$, $B_4/4=-\tfrac1{120}$, $B_6/6=\tfrac1{252}$, $B_8/8=-\tfrac1{240}$, $|B_{10}|/10=\tfrac1{132}$, and $c_4=\tfrac12-\gamma+\tfrac1{12}-\tfrac1{120}+\tfrac1{252}-\tfrac1{240}=\tfrac{2897}{5040}-\gamma$. Finally $\sum_{d\mid n,d>1}d^{-10}\le\sum_{n\ge2}n^{-10}\le2^{-10}+3^{-10}+\int_3^\infty x^{-10}dx=2^{-10}+3^{-10}+\tfrac{3^{-9}}9<10^{-3}$, so $|E_4(n)|<\tfrac1{132000}$. $\square$

*Check.* (B.2) was compared with exact rational values of $F(n)$ for all 2862 squarefree $n\le6000$ with $\omega(n)\ge2$ (and (6.1) below for all 1140 non-squarefree ones with $n\le3000$): the error never exceeded its bound (largest ratio $0.9995$, at $n=71\cdot73$, where $\theta_4(71)\approx1$ — the bound is essentially sharp).

## 5. A uniform bound and the proof of Theorem A

**Lemma 5.** For every $s>1$ and every squarefree $n$: $P_s(n)>1/\zeta(s)$. In particular $P_2>6/\pi^2$, $P_4>90/\pi^4$, $P_6>945/\pi^6$, $P_8>9450/\pi^8$.
*Proof.* $\zeta(s)=\prod_p(1-p^{-s})^{-1}$ for $s>1$ (DLMF 25.2.11), so $P_s(n)=\zeta(s)^{-1}\prod_{p\nmid n}(1-p^{-s})^{-1}>\zeta(s)^{-1}$, the product over the remaining primes lying in $(0,1)$; and $\zeta(2)=\pi^2/6$, $\zeta(4)=\pi^4/90$, $\zeta(6)=\pi^6/945$, $\zeta(8)=\pi^8/9450$ (DLMF 25.6.2). $\square$

Define
$$UB(k):=\tfrac12P_1(k)+\tfrac1{120}P_4(k)+\tfrac1{240}P_8(k)-\frac1{2\pi^2}-\frac{15}{4\pi^6}+\frac{2897}{5040}-\gamma+\frac1{132000}.$$

**Proposition 6.** (a) $G(p_k\#)<UB(k)$ for all $k\ge2$. (b) $UB(k+1)<UB(k)$ for all $k\ge1$. (c) $UB(92)<-3.94\cdot10^{-5}<0$.
*Proof.* (a) Apply (B.2) to $n=p_k\#$: by Lemma 5, $-\tfrac1{12}P_2(k)<-\tfrac1{12}\cdot\tfrac6{\pi^2}=-\tfrac1{2\pi^2}$ and $-\tfrac1{252}P_6(k)<-\tfrac1{252}\cdot\tfrac{945}{\pi^6}=-\tfrac{15}{4\pi^6}$, and $E_4<\tfrac1{132000}$. (b) Passing from $k$ to $k+1$ multiplies each of $P_1,P_4,P_8$ by a factor in $(0,1)$; the other terms are constants. (c) is the certificate below. $\square$

*Proof of Theorem A.* For $k\ge92$, Proposition 6 gives $G(p_k\#)<UB(k)\le UB(92)<0$, so by Theorem B $\operatorname{sign}a(p_k\#)=(-1)^k\operatorname{sign}G(p_k\#)=(-1)^{k+1}$. As $p_{92}=479$, $p_{93}=487$: $a(479\#)<0$ with $\omega=92$ even (so "$a(n)<0\Rightarrow\omega(n)$ odd" fails), and $a(487\#)>0$ with $\omega=93$ odd (so "$\omega(n)$ odd $\ge3\Rightarrow a(n)<0$" fails). $\square$

**Certificate for Proposition 6(c).** Only finitely many numbers enter: **(N1)** the first 92 primes $2,3,\dots,479$; **(N2)** the exact rationals $P_1(92)=\prod_{i\le92}(1-1/p_i)$, $P_4(92)$, $P_8(92)$; **(N3)** $\pi$ and $\gamma$. An upper bound for $UB(92)$ results from upper bounds for the three positive product terms, an *upper* bound $\pi\le3.1415926536$ (it makes $-1/(2\pi^2)$ and $-15/(4\pi^6)$ larger) and a *lower* bound $\gamma\ge0.5772156649$:

| term | $\tfrac12P_1(92)$ | $\tfrac1{120}P_4(92)$ | $\tfrac1{240}P_8(92)$ | $-\tfrac1{2\pi^2}$ | $-\tfrac{15}{4\pi^6}$ | $\tfrac{2897}{5040}$ | $-\gamma$ | $\tfrac1{132000}$ | **sum** |
|---|---|---|---|---|---|---|---|---|---|
| $\le$ | $0.0450790249$ | $0.0076994867$ | $0.0041497467$ | $-0.0506605918$ | $-0.0039006055$ | $0.5748015874$ | $-0.5772156649$ | $0.0000075758$ | $\mathbf{-0.0000394407}$ |

Hence $UB(92)\le-3.944\cdot10^{-5}<0$. Reference values: $P_1(92)=0.0901580497409755\ldots$, $P_4(92)=0.9239384033245916\ldots$, $P_8(92)=0.9959392011255151\ldots$, $UB(92)=-3.9440952182439\ldots\cdot10^{-5}$ (interval arithmetic, `mpmath.iv`, 30 digits: $UB(92)\in[-3.94409521824394630,-3.94409521824394629]\cdot10^{-5}$). **Precision required:** each of the eight terms to absolute accuracy $10^{-6}$ (total error $<10^{-5}$, well inside the margin $3.9\cdot10^{-5}$); the table uses 10 decimals with directed rounding, the products being exact rationals. For comparison $UB(91)=+5.49\cdot10^{-5}>0$, so the uniform bound singles out $k=92$ exactly.

## 6. Remarks

**(a) $k\le91$, and certified values.** Using the per-$n$ error bound of Theorem B instead of $1/132000$, one gets $G(p_k\#)>0$ for $2\le k\le91$ (smallest margin at $k=91$: $G(467\#)\in[2.51,4.01]\cdot10^{-5}$), and $G(479\#)\in[-6.90,-5.40]\cdot10^{-5}$, $G(487\#)\in[-1.61,-1.46]\cdot10^{-4}$. Thus $479\#$ (199 digits) is the first primorial violating the conjecture. Numbers needed: $P_s(k)$ for $s\in\{1,2,4,6,8,10\}$, $k\le93$, to $10^{-7}$.

**(b) The heuristic sign law.** Since $P_2=P_1Q$, (B.2) reads
$$G(n)=\tfrac12P_1(n)\Big(1-\frac{Q(n)}6\Big)+\Big[\tfrac1{120}P_4-\tfrac1{252}P_6+\tfrac1{240}P_8+c_4\Big]+E_4(n),$$
where the bracket lies in $[0.00546,0.00619]$ for every squarefree $n$ (use $\zeta(s)^{-1}\le P_s\le1$) and tends to $\frac1{120\zeta(4)}-\frac1{252\zeta(6)}+\frac1{240\zeta(8)}+c_4=0.0055346$ as the primes grow. So $G(n)\approx\tfrac12P_1(n)\big(1-\tfrac{Q(n)}6\big)+0.0055$: the sign flips when $Q(n)$ exceeds $6\big(1+0.011/P_1(n)\big)$, i.e. only when $n$ contains essentially all small primes. For primorials $P_1(k)Q(k)=P_2(k)\approx6/\pi^2$, giving the flip at $Q\approx6.736$; indeed $Q(91)=6.7308<6.736<Q(92)=6.7448$. For fixed $\omega$, $Q(n)$ is maximal at the primorial, so primorials are the extremal case.

**(c) Where the conjecture can fail.** *Proposition 7.* Let $\omega(n)\ge2$ and $q=n/\operatorname{rad}(n)$. If $Q(\operatorname{rad}n)\le6q$, then $G(n)>0$, i.e. $\operatorname{sign}a(n)=(-1)^{\omega(n)}$ as conjectured. Since $Q(\operatorname{rad}n)\le Q(\omega(n))$, $Q(51)=5.987<6<Q(52)=6.012$ and $Q(6480)=11.9998<12<Q(6481)$, every counterexample has $\omega(n)\ge52$, and every non-squarefree counterexample has $\omega(n)\ge6481$ (cases $\omega\le1$ are trivial: $a(1)=1$, $F(p^a)=H_{p^a}-H_{p^{a-1}}>0$).
*Proof.* Squarefree case ($q=1$): in (B.2), $\tfrac12P_1(1-Q/6)\ge0$, and with $P_4>90/\pi^4$, $P_8>9450/\pi^8$ (Lemma 5), $P_6\le1$: $G(n)\ge\frac{90}{120\pi^4}-\frac1{252}+\frac{9450}{240\pi^8}+\frac{2897}{5040}-\gamma-\frac1{132000}=0.005459\ldots>0$.
Non-squarefree case: the divisors $d$ of $n$ with $\mu(n/d)\ne0$ are $d=qe$, $e\mid\operatorname{rad}n$, with $\mu(n/d)=\mu(\operatorname{rad}(n)/e)=(-1)^{\omega(n)-\omega(e)}$; so Corollary 2, Lemma 3 (for $\operatorname{rad}n$ and $f(e)=\varepsilon(qe)$) and (3.2) give, with $P_s=P_s(\operatorname{rad}n)$,
$$G(n)=\frac{P_1}{2q}-\frac{P_2}{12q^2}+\frac{P_4}{120q^4}-\frac{P_6}{252q^6}+\frac{P_8}{240q^8}+E,\qquad |E|\le\frac{1}{132\,q^{10}}\prod_{p\mid n}(1+p^{-10})<\frac{1.001}{132\,q^{10}}\tag{6.1}$$
(no constant term: $d=1$ does not occur). If $Q\le6q$ the first two terms give $\frac{P_1}{2q}(1-\frac{Q}{6q})\ge0$, and the rest is $\ge q^{-4}\big[\frac{90}{120\pi^4}-\frac1{252q^2}-\frac{1.001}{132q^6}\big]\ge q^{-4}[0.007699-0.000992-0.000119]>0$ for $q\ge2$. $\square$

**(d) Non-squarefree counterexamples exist as well.** For $n_k:=4\cdot3\cdot5\cdots p_k=2\,p_k\#$ ($q=2$, $\operatorname{rad}n_k=p_k\#$), keep the term $e=1$ of (6.1) exactly, $R_4(2)=\varepsilon(2)-\big[\tfrac14-\tfrac1{48}+\tfrac1{1920}-\tfrac1{16128}+\tfrac1{61440}\big]=-4.6175\cdot10^{-6}$; the remaining error is below $\frac1{132}2^{-10}(\zeta(10)-1)<7.4\cdot10^{-9}$. The resulting enclosures give $G(n_k)>0$ for $2\le k\le9230$ and $G(n_{9231})<0$ ($p_{9231}=95783$), and the decreasing uniform bound $UB_2(k):=\tfrac14P_1(k)+\tfrac{P_4(k)}{1920}+\tfrac{P_8(k)}{61440}-\tfrac1{8\pi^2}-\tfrac{945}{16128\pi^6}+R_4(2)+7.4\cdot10^{-9}$ is negative for all $k\ge9231$. Hence $\operatorname{sign}a(n_k)=(-1)^{k+1}$ for all $k\ge9231$: the conjecture also fails for non-squarefree $n$, but the first such $n_k$ has $41{,}492$ digits — consistent with the OEIS b-file ($n\le411$), with Proposition 7, and with the heuristic threshold $Q\approx6q=12$.
