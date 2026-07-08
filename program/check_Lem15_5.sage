# check_Lem15_5.sage -- verification of Lemma 15.5 (lem:sin3sys) of d10runkone.tex.
#
# The transported linear system.  With bar n = ceil(n/2), un n = floor(n/2), the
# reduction produces the un n x un n system (unknowns ytil_l, l=0..un n-1)
#   sum_{l=0}^{un n-1} C(1/2, bar n + p - l) ytil_l = C(-1/2, bar n + p),  p=0..un n-1,
# and the lemma asserts its unique solution satisfies
#   sum_{l=0}^{un n-1} (-1)^l ytil_l = - C(n,2)   ( = w^T A^{-1} u ).
#
# Two exact checks over QQ:
#   (A) solve the system directly and compare sum (-1)^l ytil_l with -C(n,2);
#   (B) the explicit solution Y(x) = -sum_{j=1}^{un n} C(n,2j)(1+x)^{j-1} of the proof
#       satisfies the system and Y(-1) = -C(n,2).
#
# Run:  sage check_Lem15_5.sage N   checks n = 1..N.

import sys

def gen_binom(a, k):
    # binomial(a, k) exact; 0 for k < 0
    if k < 0:
        return QQ(0)
    return prod(QQ(a) - j for j in range(k)) / factorial(k)

Rx = PolynomialRing(QQ, 'x')
x = Rx.gen()

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("Lemma 15.5 (lem:sin3sys) of d10runkone.tex   [exact over QQ]")
    print("")
    print("(A) solve system, sum (-1)^l ytil_l  vs  -C(n,2)")
    print("(B) explicit Y(x): system residual = 0 and Y(-1) = -C(n,2)")
    print("")
    print("  n | un n | (A) sum | -C(n,2) | dif | (B) sys-ok | Y(-1)")
    print("----+------+---------+---------+-----+------------+-------")

    all_ok = True
    for n in range(1, N + 1):
        bn = (n + 1) // 2
        un = n // 2
        target = -binomial(n, 2)

        # (A) solve the linear system
        if un == 0:
            sumA = QQ(0)
        else:
            M = matrix(QQ, un, un, lambda p, l: gen_binom(QQ(1) / 2, bn + p - l))
            rhs = vector(QQ, [gen_binom(QQ(-1) / 2, bn + p) for p in range(un)])
            ytil = M.solve_right(rhs)
            sumA = sum((-1) ** l * ytil[l] for l in range(un))
        dA = sumA - target

        # (B) explicit Y(x) and its residual in the system
        Y = -sum(binomial(n, 2 * j) * (1 + x) ** (j - 1) for j in range(1, un + 1)) if un >= 1 else Rx(0)
        sys_ok = True
        prodser = None
        if un >= 1:
            # coefficients of (1+x)^{1/2} Y(x): use (1+x)^{1/2} = sum C(1/2,k) x^k
            deg_needed = n  # need [x^k] for k = bn..n-1
            half = sum(gen_binom(QQ(1) / 2, k) * x ** k for k in range(deg_needed + 1))
            Ycoeff = Y.list()
            prodser = half * Y
            for p in range(un):
                k = bn + p
                lhs = prodser[k] if k <= prodser.degree() else QQ(0)
                if lhs != gen_binom(QQ(-1) / 2, k):
                    sys_ok = False
        Ym1 = Y(-1)
        if dA != 0 or not sys_ok or Ym1 != target:
            all_ok = False
        print("%3d | %4d | %7s | %7s | %s | %10s | %s"
              % (n, un, sumA, target, dA, sys_ok, Ym1))

    print("")
    if all_ok:
        print("Resume: sum = -C(n,2) and explicit Y solves the system for n=1..%d. Formula checked." % N)
    else:
        print("Resume: some check failed. Formula NOT checked!")
