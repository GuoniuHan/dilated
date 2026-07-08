# check_Prop11_1.sage -- verification of Proposition 11.1 (prop:sec1) of d06xcos.tex.
#
# The secant-number family at s = 0.  Let f(x) = (1+x)/cos(x), so
#   a_n = n! [x^n] (1+x)/cos(x),   i.e.
#   a_{2k} = E_{2k},  a_{2k+1} = (2k+1) E_{2k}   (secant numbers).
# Then the dilated Hankel determinant HH_n(f) = det(a_{2i+j})_{0<=i,j<=n-1} is
#   HH_n(f) = 2^C(n,2) * ((n-1)!!)^2 * prod_{k=1}^{n-2} (k!!)^6.
#
# Everything is exact over QQ (no free parameter).
#
# Run:  sage check_Prop11_1.sage N
# prints, for n = 1..N, LHS(n), RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

def df(k):
    # double factorial k!!  (k>=0);  0!! = (-1)!! = empty product = 1
    return prod(range(k, 0, -2))

def build_a(M):
    # a_0..a_M from the EGF (1+x)/cos(x), exact over QQ
    Rx = PowerSeriesRing(QQ, 'x', default_prec=M + 2)
    x = Rx.gen()
    g = (1 + x) / cos(x)
    return [factorial(m) * g[m] for m in range(M + 1)]

def LHS(n, a):
    return matrix(QQ, n, n, lambda i, j: a[2 * i + j]).det()

def RHS(n):
    return (2 ** binomial(n, 2)
            * df(n - 1) ** 2
            * prod(df(k) ** 6 for k in range(1, n - 1)))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)          # max index needed is 3n-3

    print("Proposition 11.1 (prop:sec1) of d06xcos.tex   [exact over QQ]")
    print("secant numbers a = %s ..." % a[:min(3 * N + 1, 8)])
    print("")
    print("  n |               LHS(n)               |               RHS(n)               |  dif")
    print("----+------------------------------------+------------------------------------+-----")

    all_zero = True
    for n in range(1, N + 1):
        L = LHS(n, a)
        R = RHS(n)
        d = L - R
        if d != 0:
            all_zero = False
        print("%3d | %34s | %34s | %s" % (n, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
