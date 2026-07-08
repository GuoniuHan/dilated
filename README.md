# Dilated Hankel determinants

*Guo-Niu Han* — 2026/07/08


**Status.** Complete. Every displayed theorem, proposition, lemma, corollary,
and named identity of the paper (Sections 4–27 and 29) now has a standalone
SageMath checker in [`program/`](program), with a saved run in
[`output/`](output); all of them report `dif = 0`. Sections 1–3 (introduction,
preliminaries, the reduction principles) and Section 28 (the
Lindström–Gessel–Viennot discussion) state no determinant evaluation to check,
and Section 30 is the concluding remarks.

## Abstract

For a sequence $`\mathbf{a}=(a_0,a_1,\dots)`$ we define its *dilated Hankel determinant*
$`\ddot{H}_n(\mathbf{a})=\det(a_{2i+j})_{0\le i,j\le n-1}`$, the minor
of the infinite Hankel matrix $`(a_{i+j})`$ formed from the even-indexed rows and
the first $`n`$ columns. We prove that, for a broad class of sequences,
$`\ddot{H}_n`$ admits a remarkably simple product evaluation. This mirrors the
behaviour of the classical Hankel determinant $`H_n`$, but with two key
distinctions: the class of sequences for which such formulas are known is far
larger in the classical case; and, whereas $`H_n`$ enjoys a single universal
evaluation — the Heilermann formula via the Jacobi continued fraction — no
analogous general method exists for the dilated determinant, which is therefore
considerably more challenging. Our evaluations instead rest on six reduction
principles developed here, four of general scope and two of a more specialised
nature. The cases treated include the factorial numbers, the Catalan and
central binomial coefficients; the Euler numbers and a one-parameter secant
family; the involution numbers; the Springer numbers along with elliptic and
derivative deformations; the reciprocal-sine function, whose evaluation rests on
a new Catalan determinant proved by condensation; a Bessel analogue of the Euler
numbers; and a multiplicative Bessel family. As an application, we settle a
conjecture of Chapoton and the author on the roots of the Poupard and Kreweras
polynomials.

## Programs (old, keep, we will work later)

- [`r1root.sage`](r1root.sage) (Sage source) — roots of the Poupard and Kreweras
  polynomials (Chapoton–Han conjecture)
- [`r1root.txt`](r1root.txt) (output)
- [`r2coro.sage`](r2coro.sage) (Sage source) — dilated Hankel determinants for
  every corollary (Beta, Gaussian, Euler, secant, rank-one, Springer,
  reciprocal-sine, elliptic, algebraic, squared-algebraic, Bessel and
  multiplicative-Bessel families), with single and double shifts
- [`r2coro.txt`](r2coro.txt) (output)

## Check formulas

Each formula of the paper has a standalone [SageMath](https://www.sagemath.org)
program in [`program/`](program). Running

```
sage program/check_XXX.sage N
```

prints, for every index up to $`N`$ (a single $`n=1,\dots,N`$ for a determinant
evaluation, or a pair such as $`(p,q)`$ for a moment identity), the left-hand side
$`\mathrm{LHS}`$ (computed directly), the right-hand side $`\mathrm{RHS}`$ (the
claimed closed form), and their difference $`\mathrm{dif}=\mathrm{LHS}-\mathrm{RHS}`$;
it ends with a one-line résumé stating whether all differences vanish. All
arithmetic is exact (rational, or symbolic in a parameter such as $`c`$, or in a
$`\Gamma`$-value), so `dif = 0` is a genuine identity, not a numerical
approximation. The parameters ($`\rho,\alpha,\beta`$, the step $`r`$, the row indices
$`x_i`$, …) are set at the top of each file and may be freely changed. A saved run
with `N = 5` is stored alongside each program in [`output/`](output).

### Section 4 — the Beta family (`d02beta.tex`)

| Statement | Sequence / parameters | Program | Output |
|---|---|---|---|
| Theorem 4.1 (general, arbitrary nodes $`x_i`$) | $`a_m=a_0\rho^m(\alpha)_m/(\beta)_m`$ | [`check_Thm4_1.sage`](program/check_Thm4_1.sage) | [`.txt`](output/check_Thm4_1.txt) |
| Proposition 4.2 (degenerate, $`\beta`$ absent) | $`a_m=a_0\rho^m(\alpha)_m`$ | [`check_Prop4_2.sage`](program/check_Prop4_2.sage) | [`.txt`](output/check_Prop4_2.txt) |
| Corollary 4.3 (dilated Beta family) | $`a_m=a_0\rho^m(\alpha)_m/(\beta)_m`$ | [`check_Coro4_3.sage`](program/check_Coro4_3.sage) | [`.txt`](output/check_Coro4_3.txt) |
| Corollary 4.4 (factorials) | $`a_m=(m+r)!`$ | [`check_Coro4_4.sage`](program/check_Coro4_4.sage) | [`.txt`](output/check_Coro4_4.txt) |
| Corollary 4.5 (Gamma/Laguerre continuation) | $`a_m=\Gamma(m+\alpha+1)`$ | [`check_Coro4_5.sage`](program/check_Coro4_5.sage) | [`.txt`](output/check_Coro4_5.txt) |
| Corollary 4.6 ($`r`$-step Catalan) | $`a_m=C_m=\binom{2m}{m}/(m+1)`$ | [`check_Coro4_6.sage`](program/check_Coro4_6.sage) | [`.txt`](output/check_Coro4_6.txt) |
| Corollary 4.7 ($`r`$-step central binomial) | $`a_m=\binom{2m}{m}`$ | [`check_Coro4_7.sage`](program/check_Coro4_7.sage) | [`.txt`](output/check_Coro4_7.txt) |
| Identity (4.15) | $`a_m=(2m+1)!!`$ | [`check_Eq4_15.sage`](program/check_Eq4_15.sage) | [`.txt`](output/check_Eq4_15.txt) |
| Identity (4.16) | $`a_m=(2m)!/m!`$ | [`check_Eq4_16.sage`](program/check_Eq4_16.sage) | [`.txt`](output/check_Eq4_16.txt) |
| Identity (4.17) | $`a_m=\binom{2m+1}{m}`$ | [`check_Eq4_17.sage`](program/check_Eq4_17.sage) | [`.txt`](output/check_Eq4_17.txt) |

### Section 5 — the Gaussian family (`d03invo.tex`)

Here $`c`$ is kept symbolic, so every check is a polynomial identity in $`c`$. The
two lemmas are moment identities of the Gaussian functional $`\mathcal L`$, indexed
by a pair; `N` bounds both indices.

| Statement | Content | Program | Output |
|---|---|---|---|
| Theorem 5.1 (Gaussian family) | $`\ddot{H}_n=(2c)^{\binom n2}\prod_{k<n}k!`$ for $`f=e^{cx+x^2/2}`$ | [`check_Thm5_1.sage`](program/check_Thm5_1.sage) | [`.txt`](output/check_Thm5_1.txt) |
| Lemma 5.3 (mixed Hermite moments) | $`\mathcal L[\mathrm{He}_p(c+z)\mathrm{He}_q(z)]=q!\binom{p}{q}c^{p-q}`$ | [`check_Lem5_3.sage`](program/check_Lem5_3.sage) | [`.txt`](output/check_Lem5_3.txt) |
| Lemma 5.4 (triangularising family) | $`\mathcal L[\Phi_i(z)\mathrm{He}_j(z)]=j!\binom{i}{j-i}(2c)^{2i-j}`$ | [`check_Lem5_4.sage`](program/check_Lem5_4.sage) | [`.txt`](output/check_Lem5_4.txt) |

### Section 6 — the single shift of the Gaussian family (`d03invo.tex`)

The single shift is the dilated Hankel determinant of $`f'`$, i.e.
$`\det(a_{2i+j+1})`$. Again $`c`$ is symbolic; identity (6.1) is a moment
identity, checked over the triangular set $`0\le j\le i\le N`$ (the cases it
claims).

| Statement | Content | Program | Output |
|---|---|---|---|
| Theorem 6.1 (single shift) | $`\ddot{H}_n(f')=c^n(2c)^{\binom n2}\prod_{k<n}k!`$ | [`check_Thm6_1.sage`](program/check_Thm6_1.sage) | [`.txt`](output/check_Thm6_1.txt) |
| Identity (6.1) | $`\mathcal L[\Psi_i\,u\,\mathrm{He}_j]=0\ (j<i),\ c\,i!\,(2c)^i\ (j=i)`$ | [`check_Eq6_1.sage`](program/check_Eq6_1.sage) | [`.txt`](output/check_Eq6_1.txt) |

### Section 7 — the Euler number family (`d04euler.tex`)

| Statement | Program | Output |
|---|---|---|
| Proposition 7.1 (`prop:family`) | [`check_Prop7_1.sage`](program/check_Prop7_1.sage) | [`.txt`](output/check_Prop7_1.txt) |
| Lemma 7.2 (`lem:cf`) | [`check_Lem7_2.sage`](program/check_Lem7_2.sage) | [`.txt`](output/check_Lem7_2.txt) |
| Lemma 7.3 (`lem:connGene`) | [`check_Lem7_3.sage`](program/check_Lem7_3.sage) | [`.txt`](output/check_Lem7_3.txt) |
| Lemma 7.4 (`lem:bindetSigma`) | [`check_Lem7_4.sage`](program/check_Lem7_4.sage) | [`.txt`](output/check_Lem7_4.txt) |
| Lemma 7.5 (`lem:Pahalf`) | [`check_Lem7_5.sage`](program/check_Lem7_5.sage) | [`.txt`](output/check_Lem7_5.txt) |

### Section 8 — the double shift of the Euler number family (`d04euler.tex`)

The double shift is $`\ddot H_n^{(2)}=\det(a_{2i+j+2})`$ of the Euler $`(s,t)`$
family, whose sequence $`a`$ (with $`a_{2k}=E^{(s+1)}_{2k}`$, $`a_{2k+1}=E^{(t+1)}_{2k}`$)
is generated from the Stieltjes fractions $`u_j=j(j+s)`$, $`v_j=j(j+t)`$. Here
$`\bar n=\lceil n/2\rceil`$, $`\underline n=\lfloor n/2\rfloor`$, and $`\ddot H_n`$
is the unshifted determinant of the same $`a`$. Parameters $`s,t`$ are set at the
top; keep $`(t-s)/2`$ away from the non-negative integers or the family degenerates.

| Statement | Content | Program | Output |
|---|---|---|---|
| Proposition 8.4 (double shift, special family) | $`\ddot H_n^{(2)}=(2n-1)!!\prod_{k=0}^{\bar n-1}(s+2k+1)\prod_{k=0}^{\underline n-1}(t+2k+1)\,\ddot H_n`$ | [`check_Prop8_4.sage`](program/check_Prop8_4.sage) | [`.txt`](output/check_Prop8_4.txt) |

### Section 9 — the single shift of the Euler number family on the line `t=1` (`d05eulershift.tex`)

| Statement | Program | Output |
|---|---|---|
| Proposition 9.1 (`prop:shift`) | [`check_Prop9_1.sage`](program/check_Prop9_1.sage) | [`.txt`](output/check_Prop9_1.txt) |
| Lemma 9.2 (`lem:mixedconn`) | [`check_Lem9_2.sage`](program/check_Lem9_2.sage) | [`.txt`](output/check_Lem9_2.txt) |
| equation (`9.3`) (`eq:shiftmonic`) | [`check_Eq9_3.sage`](program/check_Eq9_3.sage) | [`.txt`](output/check_Eq9_3.txt) |
| equation (`9.4`) (`eq:shiftconst`) | [`check_Eq9_4.sage`](program/check_Eq9_4.sage) | [`.txt`](output/check_Eq9_4.txt) |

### Section 10 — the single shift of the Euler family on the line $`t=3`$ (`d05eulershift.tex`)

Here $`t=3`$ and $`s`$ is symbolic. The single shift is
$`\ddot H_n^{(1)}=\det(a_{2i+j+1})`$, built from the moments of the Stieltjes
fractions $`u_j=j(j+s)`$ (even) and $`v_j=j(j+3)`$ (odd); the orthogonal
polynomials are generated from their three-term recurrences. Write
$`\bar n=\lceil n/2\rceil`$, $`\underline n=\lfloor n/2\rfloor`$,
$`b=(s-3)/2`$, and $`\tilde\kappa`$ for the mixed connection coefficients.
Every determinant identity is checked as a polynomial identity in $`s`$; the
divergence/rank statements are checked at the relevant $`s=-(2k+1)`$.

| Statement | Content | Program | Output |
|---|---|---|---|
| Lemma 10.1 (mixed connection at $`t=3`$) | $`\tilde\kappa_{i,m}=\tfrac{(2i)!}{(2m)!}\binom{(s-1)/2}{i-m}+(s+1)\tfrac{(2i)!}{(2m+1)!}\binom{(s-3)/2}{i-m-1}`$ | [`check_Lem10_1.sage`](program/check_Lem10_1.sage) | [`.txt`](output/check_Lem10_1.txt) |
| Lemma 10.3 (Cauchy–Binet collapse) | $`\det(\tilde\kappa_{\bar n+r,m})=(\text{factor})\det N`$, $`N=VB`$, $`\det N=\sum_p c_p D_p`$ | [`check_Lem10_3.sage`](program/check_Lem10_3.sage) | [`.txt`](output/check_Lem10_3.txt) |
| Proposition 10.5 (closed evaluation) | $`\ddot H_n^{(1)}=(-1)^{\binom{\bar n}2}(\prod h^T_l)(\prod h^{S'}_c/(2c+1)!)(\prod(2\bar n+2r)!)\sum_p c_p D_p`$ | [`check_Prop10_5.sage`](program/check_Prop10_5.sage) | [`.txt`](output/check_Prop10_5.txt) |
| Theorem 10.6 (factorisation; carrier) | $`\ddot H_n^{(1)}`$ splits into linear forms over $`\mathbb Q`$; unique non-integer factor $`\Gamma_n=(2\bar n+1)s+(4\bar n\underline n-2\bar n-1)`$ | [`check_Thm10_6.sage`](program/check_Thm10_6.sage) | [`.txt`](output/check_Thm10_6.txt) |
| Lemma 10.7 (atomic degeneration) | $`a_{2j}\rvert_{s=-(2k+1)}=\sum_l\omega_l(-4l^2)^j`$; $`(s+2k+1)^{\underline n-k}\mid\ddot H_n^{(1)}`$ | [`check_Lem10_7.sage`](program/check_Lem10_7.sage) | [`.txt`](output/check_Lem10_7.txt) |
| Lemma 10.8 (rank of $`N`$) | $`\operatorname{rank}N\rvert_{s=-(2k+1)}\le k+1`$; $`(s+2k+1)^{\underline n-1-k}\mid\det N`$ | [`check_Lem10_8.sage`](program/check_Lem10_8.sage) | [`.txt`](output/check_Lem10_8.txt) |
| Identity (10.8) (Chu–Vandermonde) | $`R=\sum_{j=0}^{k+1}\binom{k+X}{j}G_j(m)`$, and $`N_{r,m}\rvert_{s=-(2k+1)}=(-1)^uR`$ | [`check_Eq10_8.sage`](program/check_Eq10_8.sage) | [`.txt`](output/check_Eq10_8.txt) |
| Lemma 10.10 (trace identity) | $`\operatorname{tr}(\bar M^{-1}D_a\bar M D_m)=\sum_r a_r(\bar n+r)-\bar n a_0\underline n`$ | [`check_Lem10_10.sage`](program/check_Lem10_10.sage) | [`.txt`](output/check_Lem10_10.txt) |
| Lemma 10.11 (contiguous relation) | $`a_{2k}(s)=\tfrac{a_{2k+2}(s-2)+(s-1)^2a_{2k}(s-2)}{s(s-1)}`$; $`(s(s-1))^{\underline n}\ddot H_n^{(1)}(s)=\sum_q(s-1)^{2q}G_q(s-2)`$ | [`check_Lem10_11.sage`](program/check_Lem10_11.sage) | [`.txt`](output/check_Lem10_11.txt) |

### Section 11 — the secant-number family `(1+x)/cos^{s+1}x`, case `s=0` (`d06xcos.tex`)

| Statement | Program | Output |
|---|---|---|
| Proposition 11.1 (`prop:sec1`) | [`check_Prop11_1.sage`](program/check_Prop11_1.sage) | [`.txt`](output/check_Prop11_1.txt) |
| Lemma 11.2 (`lem:cdh`) | [`check_Lem11_2.sage`](program/check_Lem11_2.sage) | [`.txt`](output/check_Lem11_2.txt) |
| Lemma 11.4 (`lem:conn3`) | [`check_Lem11_4.sage`](program/check_Lem11_4.sage) | [`.txt`](output/check_Lem11_4.txt) |

### Section 12 — the secant-number family `(1+x)/cos^{s+1}x`, general `s` (`d06xcos.tex`)

| Statement | Program | Output |
|---|---|---|
| Theorem 12.1 (`thm:allstar`) | [`check_Thm12_1.sage`](program/check_Thm12_1.sage) | [`.txt`](output/check_Thm12_1.txt) |
| Lemma 12.2 (`lem:reduction`) | [`check_Lem12_2.sage`](program/check_Lem12_2.sage) | [`.txt`](output/check_Lem12_2.txt) |
| Lemma 12.3 (`lem:rankdrop`) | [`check_Lem12_3.sage`](program/check_Lem12_3.sage) | [`.txt`](output/check_Lem12_3.txt) |
| Lemma 12.4 (`lem:deriv`) | [`check_Lem12_4.sage`](program/check_Lem12_4.sage) | [`.txt`](output/check_Lem12_4.txt) |

### Section 13 — the single shift of the secant-number family (`d07xcosshift.tex`)

| Statement | Program | Output |
|---|---|---|
| Theorem 13.1 (`thm:shift`) | [`check_Thm13_1.sage`](program/check_Thm13_1.sage) | [`.txt`](output/check_Thm13_1.txt) |
| Lemma 13.3 (`lem:oddred`) | [`check_Lem13_3.sage`](program/check_Lem13_3.sage) | [`.txt`](output/check_Lem13_3.txt) |
| Lemma 13.4 (`lem:div`) | [`check_Lem13_4.sage`](program/check_Lem13_4.sage) | [`.txt`](output/check_Lem13_4.txt) |
| Lemma 13.5 (`lem:percol`) | [`check_Lem13_5.sage`](program/check_Lem13_5.sage) | [`.txt`](output/check_Lem13_5.txt) |
| Lemma 13.6 (`lem:deg`) | [`check_Lem13_6.sage`](program/check_Lem13_6.sage) | [`.txt`](output/check_Lem13_6.txt) |
| Proposition 13.7 (`prop:const`) | [`check_Prop13_7.sage`](program/check_Prop13_7.sage) | [`.txt`](output/check_Prop13_7.txt) |
| Lemma 13.8 (`lem:samekappa`) | [`check_Lem13_8.sage`](program/check_Lem13_8.sage) | [`.txt`](output/check_Lem13_8.txt) |

### Section 14 — the double shift of the secant-number family at `s=1` (`d08xcosds.tex`)

| Statement | Program | Output |
|---|---|---|
| Theorem 14.1 (`thm:dblshift`) | [`check_Thm14_1.sage`](program/check_Thm14_1.sage) | [`.txt`](output/check_Thm14_1.txt) |
| Lemma 14.2 (`lem:ds-red`) | [`check_Lem14_2.sage`](program/check_Lem14_2.sage) | [`.txt`](output/check_Lem14_2.txt) |
| Lemma 14.3 (`lem:ds-structred`) | [`check_Lem14_3.sage`](program/check_Lem14_3.sage) | [`.txt`](output/check_Lem14_3.txt) |
| Lemma 14.4 (`lem:ds-Delta`) | [`check_Lem14_4.sage`](program/check_Lem14_4.sage) | [`.txt`](output/check_Lem14_4.txt) |
| Lemma 14.5 (`lem:ds-cauchy`) | [`check_Lem14_5.sage`](program/check_Lem14_5.sage) | [`.txt`](output/check_Lem14_5.txt) |

### Section 15 — a rank-one perturbation of the Euler number family (`d10runkone.tex`)

| Statement | Program | Output |
|---|---|---|
| Proposition 15.1 (`prop:sin3`) | [`check_Prop15_1.sage`](program/check_Prop15_1.sage) | [`.txt`](output/check_Prop15_1.txt) |
| Lemma 15.2 (`lem:rank1`) | [`check_Lem15_2.sage`](program/check_Lem15_2.sage) | [`.txt`](output/check_Lem15_2.txt) |
| Lemma 15.3 (`lem:cfshift`) | [`check_Lem15_3.sage`](program/check_Lem15_3.sage) | [`.txt`](output/check_Lem15_3.txt) |
| Lemma 15.4 (`lem:evalminus1`) | [`check_Lem15_4.sage`](program/check_Lem15_4.sage) | [`.txt`](output/check_Lem15_4.txt) |
| Lemma 15.5 (`lem:sin3sys`) | [`check_Lem15_5.sage`](program/check_Lem15_5.sage) | [`.txt`](output/check_Lem15_5.txt) |

### Section 16 — the Springer number family (`d11springer.tex`)

The Springer family is $`1/(\cos x-t\sin x)^r`$; its moments are
$`a_n(t)=n!\,[x^n](\cos x-t\sin x)^{-r}`$. The exponent $`r\ge1`$ is a parameter
at the top; $`t`$ is symbolic, so the identity is checked in $`\mathbb Q[t]`$ and
the LHS is printed factored, exhibiting the $`(t(t^2+1))^{\binom n2}`$ divisor.
($`r=t=1`$ is the Springer case $`\ddot H_n=4^{\binom n2}\prod k!(2k)!`$.)

| Statement | Content | Program | Output |
|---|---|---|---|
| Theorem 16.1 (Springer family) | $`\ddot H_n\!\big(\tfrac{1}{(\cos x-t\sin x)^r}\big)=(t(t^2+1))^{\binom n2}\,\ddot H_n\!\big(\tfrac{1}{(1-x)^r}\big)`$ | [`check_Thm16_1.sage`](program/check_Thm16_1.sage) | [`.txt`](output/check_Thm16_1.txt) |

### Section 17 — a derivative of the Springer number family at `t=1` (`d11springer.tex`)

| Statement | Program | Output |
|---|---|---|
| Theorem 17.1 (`thm:deriv`) | [`check_Thm17_1.sage`](program/check_Thm17_1.sage) | [`.txt`](output/check_Thm17_1.txt) |
| Lemma 17.2 (`lem:Apoly`) | [`check_Lem17_2.sage`](program/check_Lem17_2.sage) | [`.txt`](output/check_Lem17_2.txt) |
| Lemma 17.3 (`lem:gausskernel`) | [`check_Lem17_3.sage`](program/check_Lem17_3.sage) | [`.txt`](output/check_Lem17_3.txt) |

### Section 18 — the reciprocal-sine case $`(1+x)\,x/\sin x`$ (`d12xsinx.tex`)

With $`g=x/\sin x`$, $`b_k=(2k)![x^{2k}]g`$, the moments are $`a_{2k}=b_k`$,
$`a_{2k+1}=(2k+1)b_k`$. The proof runs through the connection coefficients
$`\kappa_{i,m}`$ between the two Wilson functionals and a Catalan determinant.
All checks use exact rational arithmetic.

| Statement | Content | Program | Output |
|---|---|---|---|
| Theorem 18.1 (closed form) | $`\ddot H_{n+2}=\tfrac23 Q(n)\ddot H_n`$, $`\ddot H_{2N}=(\tfrac23)^N\prod_{k<N}Q(2k)`$, $`\ddot H_{2N+1}=(\tfrac23)^N\prod_{k<N}Q(2k+1)`$ | [`check_Thm18_1.sage`](program/check_Thm18_1.sage) | [`.txt`](output/check_Thm18_1.txt) |
| Lemma 18.6 (closed form of $`\kappa`$) | $`\kappa_{i,m}=\tfrac{(3/2-d)_d(i+m+2)_d((m+1)_d)^3}{d!\,(m+3/2)_d(i+m+1/2)_d}`$, $`d=i-m`$ | [`check_Lem18_6.sage`](program/check_Lem18_6.sage) | [`.txt`](output/check_Lem18_6.txt) |
| Lemma 18.7 (Catalan splitting) | $`\kappa_{i,m}=\alpha_i\beta_m C_{i-m-1}C_{i+m}`$ (for $`i>m`$) | [`check_Lem18_7.sage`](program/check_Lem18_7.sage) | [`.txt`](output/check_Lem18_7.txt) |
| Theorem 18.8 (Catalan determinant) | $`F_N(a,b)=\det(C_{a+r-m}C_{b+r+m})`$ equals the staircase product (eq:xsin-catalan), for $`a\ge N-1,\ b\ge a+1`$ | [`check_Thm18_8.sage`](program/check_Thm18_8.sage) | [`.txt`](output/check_Thm18_8.txt) |

### Section 19 — an elliptic deformation of the Euler numbers (`d15elliptic.tex`)

| Statement | Program | Output |
|---|---|---|
| Theorem 19.1 (`thm:ell`) | [`check_Thm19_1.sage`](program/check_Thm19_1.sage) | [`.txt`](output/check_Thm19_1.txt) |

### Section 20 — an algebraic family `(1+x)/(1-x^2)^{s/2}` (`d18alge.tex`)

| Statement | Program | Output |
|---|---|---|
| Theorem 20.1 (`conj:alg`) | [`check_Thm20_1.sage`](program/check_Thm20_1.sage) | [`.txt`](output/check_Thm20_1.txt) |

### Section 21 — the squared algebraic family `(1+x)^2/(1-x^2)^{s/2}` (`d19algesq.tex`)

| Statement | Program | Output |
|---|---|---|
| Theorem 21.1 (`conj:algsq`) | [`check_Thm21_1.sage`](program/check_Thm21_1.sage) | [`.txt`](output/check_Thm21_1.txt) |

### Section 22 — shifted determinants of the two algebraic families (`d19algesq.tex`)

| Statement | Program | Output |
|---|---|---|
| Proposition 22.1 (`prop:algshift`) | [`check_Prop22_1.sage`](program/check_Prop22_1.sage) | [`.txt`](output/check_Prop22_1.txt) |
| Theorem 22.3 (`thm:algsqshift`) | [`check_Thm22_3.sage`](program/check_Thm22_3.sage) | [`.txt`](output/check_Thm22_3.txt) |

### Section 23 — a Bessel analogue of the Euler number family (`d24Bessel.tex`)

| Statement | Program | Output |
|---|---|---|
| Theorem 23.1 (`thm:besselst`) | [`check_Thm23_1.sage`](program/check_Thm23_1.sage) | [`.txt`](output/check_Thm23_1.txt) |
| Example 23.3 (`ex:besselJ0Si`) | [`check_Ex23_3.sage`](program/check_Ex23_3.sage) | [`.txt`](output/check_Ex23_3.txt) |
| Lemma 23.6 (`lem:bessel-kappa`) | [`check_Lem23_6.sage`](program/check_Lem23_6.sage) | [`.txt`](output/check_Lem23_6.txt) |
| Theorem 23.7 (`thm:bessel-det`) | [`check_Thm23_7.sage`](program/check_Thm23_7.sage) | [`.txt`](output/check_Thm23_7.txt) |

### Section 24 — the double shift of the Bessel `(s,t)` family (`d24Bessel.tex`)

| Statement | Program | Output |
|---|---|---|
| Proposition 24.2 (`prop:besseldshift`) | [`check_Prop24_2.sage`](program/check_Prop24_2.sage) | [`.txt`](output/check_Prop24_2.txt) |

### Section 25 — a multiplicative Bessel family `(1+x)cosb_v^2` at even orders (`d25xBessel.tex`)

| Statement | Program | Output |
|---|---|---|
| Lemma 25.1 (`lem:xbe-moments`) | [`check_Lem25_1.sage`](program/check_Lem25_1.sage) | [`.txt`](output/check_Lem25_1.txt) |
| Proposition 25.2 (`prop:xbe-vdm`) | [`check_Prop25_2.sage`](program/check_Prop25_2.sage) | [`.txt`](output/check_Prop25_2.txt) |
| Lemma 25.3 (`lem:xbe-coeffdet`) | [`check_Lem25_3.sage`](program/check_Lem25_3.sage) | [`.txt`](output/check_Lem25_3.txt) |
| Theorem 25.4 (`thm:xbesseleven`) | [`check_Thm25_4.sage`](program/check_Thm25_4.sage) | [`.txt`](output/check_Thm25_4.txt) |
| Corollary 25.6 (`cor:xbe-sinc`) | [`check_Coro25_6.sage`](program/check_Coro25_6.sage) | [`.txt`](output/check_Coro25_6.txt) |

### Section 26 — connection with a determinant of Chapoton--Han (`d50chapoton.tex`)

| Statement | Program | Output |
|---|---|---|
| Theorem 26.1 (`thm:rootlink`) | [`check_Thm26_1.sage`](program/check_Thm26_1.sage) | [`.txt`](output/check_Thm26_1.txt) |
| Lemma 26.2 (`lem:rhored`) | [`check_Lem26_2.sage`](program/check_Lem26_2.sage) | [`.txt`](output/check_Lem26_2.txt) |
| Theorem 26.3 (`thm:conj54double`) | [`check_Thm26_3.sage`](program/check_Thm26_3.sage) | [`.txt`](output/check_Thm26_3.txt) |
| Lemma 26.5 (`lem:rhomom`) | [`check_Lem26_5.sage`](program/check_Lem26_5.sage) | [`.txt`](output/check_Lem26_5.txt) |

### Section 27 — Corollaries: explicit evaluations (`d51coro.tex`)

The compendium of explicit members. Each checker computes the sequence from its generating function (or the stated closed form of its moments), forms `det(a_{2i+j})` (and the single/double shifts where the corollary records them), and compares to the stated product; all use exact arithmetic. (The reciprocal-sine subsection has no corollary of its own — see Section 18.)

**The Beta family (`d02beta.tex`, via Cor 4.3 / Prop 4.2).**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.1 | Catalan numbers $`C_m`$ | [`check_Coro27_1.sage`](program/check_Coro27_1.sage) | [`.txt`](output/check_Coro27_1.txt) |
| 27.2 | central binomial $`\binom{2m}{m}`$ | [`check_Coro27_2.sage`](program/check_Coro27_2.sage) | [`.txt`](output/check_Coro27_2.txt) |
| 27.3 | $`\binom{2m+1}{m}`$ | [`check_Coro27_3.sage`](program/check_Coro27_3.sage) | [`.txt`](output/check_Coro27_3.txt) |
| 27.4 | factorials $`m!`$ | [`check_Coro27_4.sage`](program/check_Coro27_4.sage) | [`.txt`](output/check_Coro27_4.txt) |
| 27.5 | double factorials $`(2m\pm1)!!`$ | [`check_Coro27_5.sage`](program/check_Coro27_5.sage) | [`.txt`](output/check_Coro27_5.txt) |
| 27.6 | $`(2m)!/m!`$ | [`check_Coro27_6.sage`](program/check_Coro27_6.sage) | [`.txt`](output/check_Coro27_6.txt) |

**The Gaussian family (`d03invo.tex`).**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.7 | $`c=1`$ (involution numbers) | [`check_Coro27_7.sage`](program/check_Coro27_7.sage) | [`.txt`](output/check_Coro27_7.txt) |
| 27.8 | $`c=-1`$ (signed involutions) | [`check_Coro27_8.sage`](program/check_Coro27_8.sage) | [`.txt`](output/check_Coro27_8.txt) |

**The Euler number family (`d04euler.tex`), by $`(s,t)`$.**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.9 | $`(0,1)`$ Euler numbers | [`check_Coro27_9.sage`](program/check_Coro27_9.sage) | [`.txt`](output/check_Coro27_9.txt) |
| 27.10 | $`(1,2)`$ | [`check_Coro27_10.sage`](program/check_Coro27_10.sage) | [`.txt`](output/check_Coro27_10.txt) |
| 27.11 | $`(2,3)`$ | [`check_Coro27_11.sage`](program/check_Coro27_11.sage) | [`.txt`](output/check_Coro27_11.txt) |
| 27.12 | $`(1,0)`$ | [`check_Coro27_12.sage`](program/check_Coro27_12.sage) | [`.txt`](output/check_Coro27_12.txt) |
| 27.13 | $`(2,1)`$ | [`check_Coro27_13.sage`](program/check_Coro27_13.sage) | [`.txt`](output/check_Coro27_13.txt) |
| 27.14 | $`(3,2)`$ | [`check_Coro27_14.sage`](program/check_Coro27_14.sage) | [`.txt`](output/check_Coro27_14.txt) |
| 27.15 | $`(0,3)`$ | [`check_Coro27_15.sage`](program/check_Coro27_15.sage) | [`.txt`](output/check_Coro27_15.txt) |
| 27.16 | $`(3,0)`$ | [`check_Coro27_16.sage`](program/check_Coro27_16.sage) | [`.txt`](output/check_Coro27_16.txt) |

**The secant family $`(1+x)/\cos^{s+1}x`$ (`d06xcos.tex`).**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.17 | $`s=0`$ | [`check_Coro27_17.sage`](program/check_Coro27_17.sage) | [`.txt`](output/check_Coro27_17.txt) |
| 27.18 | $`s=1`$ | [`check_Coro27_18.sage`](program/check_Coro27_18.sage) | [`.txt`](output/check_Coro27_18.txt) |
| 27.19 | $`s=2`$ | [`check_Coro27_19.sage`](program/check_Coro27_19.sage) | [`.txt`](output/check_Coro27_19.txt) |
| 27.20 | $`s=3`$ | [`check_Coro27_20.sage`](program/check_Coro27_20.sage) | [`.txt`](output/check_Coro27_20.txt) |

**Rank-one perturbation $`(\sin x+1)/\cos^2x+s\sin x`$ (`d10runkone.tex`).**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.21 | $`s=-2`$ | [`check_Coro27_21.sage`](program/check_Coro27_21.sage) | [`.txt`](output/check_Coro27_21.txt) |
| 27.22 | $`s=-1`$ | [`check_Coro27_22.sage`](program/check_Coro27_22.sage) | [`.txt`](output/check_Coro27_22.txt) |
| 27.23 | $`s=0`$ | [`check_Coro27_23.sage`](program/check_Coro27_23.sage) | [`.txt`](output/check_Coro27_23.txt) |
| 27.24 | $`s=1`$ | [`check_Coro27_24.sage`](program/check_Coro27_24.sage) | [`.txt`](output/check_Coro27_24.txt) |
| 27.25 | $`s=2`$ | [`check_Coro27_25.sage`](program/check_Coro27_25.sage) | [`.txt`](output/check_Coro27_25.txt) |

**The Springer number family $`1/(\cos x-t\sin x)^r`$ (`d11springer.tex`).**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.26 | $`t=1,r=1`$ Springer numbers | [`check_Coro27_26.sage`](program/check_Coro27_26.sage) | [`.txt`](output/check_Coro27_26.txt) |
| 27.27 | $`t=1,r=2`$ | [`check_Coro27_27.sage`](program/check_Coro27_27.sage) | [`.txt`](output/check_Coro27_27.txt) |

**Derivative of the Springer family $`(\cos x+\sin x)/(\cos x-\sin x)^s`$ (`d11springer.tex`).**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.28 | $`s=1`$ | [`check_Coro27_28.sage`](program/check_Coro27_28.sage) | [`.txt`](output/check_Coro27_28.txt) |
| 27.29 | $`s=2`$ | [`check_Coro27_29.sage`](program/check_Coro27_29.sage) | [`.txt`](output/check_Coro27_29.txt) |
| 27.30 | $`s=3`$ | [`check_Coro27_30.sage`](program/check_Coro27_30.sage) | [`.txt`](output/check_Coro27_30.txt) |

**Elliptic deformation $`(1+\mathrm{sn})/\mathrm{cn}`$ (`d15elliptic.tex`).**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.31 | $`m=-1`$ | [`check_Coro27_31.sage`](program/check_Coro27_31.sage) | [`.txt`](output/check_Coro27_31.txt) |
| 27.32 | $`m=2`$ | [`check_Coro27_32.sage`](program/check_Coro27_32.sage) | [`.txt`](output/check_Coro27_32.txt) |
| 27.33 | $`m=1/2`$ | [`check_Coro27_33.sage`](program/check_Coro27_33.sage) | [`.txt`](output/check_Coro27_33.txt) |

**Algebraic family $`(1+x)/(1-x^2)^{s/2}`$ (`d18alge.tex`).**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.34 | $`s=-3`$ | [`check_Coro27_34.sage`](program/check_Coro27_34.sage) | [`.txt`](output/check_Coro27_34.txt) |
| 27.35 | $`s=-1`$ | [`check_Coro27_35.sage`](program/check_Coro27_35.sage) | [`.txt`](output/check_Coro27_35.txt) |
| 27.36 | $`s=1`$ | [`check_Coro27_36.sage`](program/check_Coro27_36.sage) | [`.txt`](output/check_Coro27_36.txt) |
| 27.37 | $`s=2`$ | [`check_Coro27_37.sage`](program/check_Coro27_37.sage) | [`.txt`](output/check_Coro27_37.txt) |
| 27.38 | $`s=3`$ | [`check_Coro27_38.sage`](program/check_Coro27_38.sage) | [`.txt`](output/check_Coro27_38.txt) |

**Squared-algebraic family $`(1+x)^2/(1-x^2)^{s/2}`$ (`d19algesq.tex`).**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.39 | $`s=2`$ | [`check_Coro27_39.sage`](program/check_Coro27_39.sage) | [`.txt`](output/check_Coro27_39.txt) |
| 27.40 | $`s=4`$ | [`check_Coro27_40.sage`](program/check_Coro27_40.sage) | [`.txt`](output/check_Coro27_40.txt) |

**The Bessel $`(s,t)`$ family (`d24Bessel.tex`).**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.41 | $`(0,1/2)`$, $`J_0+\mathrm{Si}`$ | [`check_Coro27_41.sage`](program/check_Coro27_41.sage) | [`.txt`](output/check_Coro27_41.txt) |
| 27.42 | $`(1/2,0)`$ | [`check_Coro27_42.sage`](program/check_Coro27_42.sage) | [`.txt`](output/check_Coro27_42.txt) |

**Multiplicative Bessel family $`(1+x)\,\mathrm{cosb}_\nu^2`$ (`d25xBessel.tex`), even orders.**

| Corollary | Member | Program | Output |
|---|---|---|---|
| 27.43 | $`\nu=0`$ | [`check_Coro27_43.sage`](program/check_Coro27_43.sage) | [`.txt`](output/check_Coro27_43.txt) |
| 27.44 | $`\nu=1/2`$ | [`check_Coro27_44.sage`](program/check_Coro27_44.sage) | [`.txt`](output/check_Coro27_44.txt) |
| 27.45 | $`\nu=1`$ | [`check_Coro27_45.sage`](program/check_Coro27_45.sage) | [`.txt`](output/check_Coro27_45.txt) |

### Section 29 — classical Hankel determinants (`d57hankel.tex`)

| Statement | Program | Output |
|---|---|---|
| Proposition 29.1 (`prop:hankel-secodd`) | [`check_Prop29_1.sage`](program/check_Prop29_1.sage) | [`.txt`](output/check_Prop29_1.txt) |
| Proposition 29.2 (`prop:hankel-sineodd`) | [`check_Prop29_2.sage`](program/check_Prop29_2.sage) | [`.txt`](output/check_Prop29_2.txt) |
| Proposition 29.3 (`prop:hankel-sineeven`) | [`check_Prop29_3.sage`](program/check_Prop29_3.sage) | [`.txt`](output/check_Prop29_3.txt) |
| Proposition 29.4 (`prop:hankel-bessel`) | [`check_Prop29_4.sage`](program/check_Prop29_4.sage) | [`.txt`](output/check_Prop29_4.txt) |
| Proposition 29.5 (`prop:hankel-springer`) | [`check_Prop29_5.sage`](program/check_Prop29_5.sage) | [`.txt`](output/check_Prop29_5.txt) |
