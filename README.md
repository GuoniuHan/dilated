# Dilated Hankel determinants

*Guo-Niu Han* — 2026/07/08


**Status.** The intention is to provide programs that check every formula
displayed in the paper. As of today Section 4 (the Beta family) is complete and
Section 5 (the Gaussian family) is under way; the remaining sections are still to
be done.

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

