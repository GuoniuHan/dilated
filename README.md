# Dilated Hankel determinants

*Guo-Niu Han*

**Status.** The intention is to provide programs that check every formula
displayed in the paper. As of today this is only partially done; the remaining
formulas are still to be completed.

## Abstract

For a sequence $\mathbf{a}=(a_0,a_1,\dots)$ we define its *dilated Hankel
determinant* $\ddot H_n(\mathbf{a})=\det(a_{2i+j})_{0\le i,j\le n-1}$, the minor
of the infinite Hankel matrix $(a_{i+j})$ formed from the even-indexed rows and
the first $n$ columns. We prove that, for a broad class of sequences,
$\ddot H_n$ admits a remarkably simple product evaluation. This mirrors the
behaviour of the classical Hankel determinant $H_n$, but with two key
distinctions: the class of sequences for which such formulas are known is far
larger in the classical case; and, whereas $H_n$ enjoys a single universal
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

## Programs

- [`r1root.sage`](r1root.sage) (Sage source) — roots of the Poupard and Kreweras
  polynomials (Chapoton–Han conjecture)
- [`r1root.txt`](r1root.txt) (output)
- [`r2coro.sage`](r2coro.sage) (Sage source) — dilated Hankel determinants for
  every corollary (Beta, Gaussian, Euler, secant, rank-one, Springer,
  reciprocal-sine, elliptic, algebraic, squared-algebraic, Bessel and
  multiplicative-Bessel families), with single and double shifts
- [`r2coro.txt`](r2coro.txt) (output)
