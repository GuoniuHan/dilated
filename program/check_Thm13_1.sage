# check_Thm13_1.sage -- verification of Theorem 13.1 (thm:shift) of d07xcosshift.tex.
#
# Single shift of the secant-number family f(x) = (1+x)/cos(x)^{s+1}.  With
#   a_{2k} = E^{(s+1)}_{2k},  a_{2k+1} = (2k+1) E^{(s+1)}_{2k}
# (equivalently a_m = m! [x^m] (1+x) cos(x)^{-(s+1)}), the single shift is
#   HH_n^{(1)} = det(a_{2i+j+1})_{0<=i,j<n},   HH_n(f) = det(a_{2i+j}).
# The theorem asserts
#   HH_n^{(1)} = ((n-1)!!)^2 * HH_n(f).
#
# The parameter s is kept SYMBOLIC over QQ, so each identity is a polynomial
# identity in s.  LHS(n) = HH_n^{(1)}; RHS(n) = ((n-1)!!)^2 HH_n(f).
#
# Run:  sage check_Thm13_1.sage N
# prints, for n = 1..N, dif = LHS(n)-RHS(n), then a resume.

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()

def df(k):
    return prod(range(k, 0, -2))          # double factorial k!!

def build_a(M):
    # a_0..a_M as polynomials in s, from the EGF (1+x) cos(x)^{-(s+1)}
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2)
    x = Rx.gen()
    g = (1 + x) * (-(s + 1) * cos(x).log()).exp()
    return [factorial(m) * g[m] for m in range(M + 1)]

def HH_shift(n, a):
    return matrix(Rs, n, n, lambda i, j: a[2 * i + j + 1]).det()

def HH(n, a):
    return matrix(Rs, n, n, lambda i, j: a[2 * i + j]).det()

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)          # max index needed is 2(n-1)+(n-1)+1 = 3n-2

    print("Theorem 13.1 (thm:shift) of d07xcosshift.tex   [s kept symbolic]")
    print("HH_n^(1) = ((n-1)!!)^2 HH_n(f)")
    print("")
    print("  n | ((n-1)!!)^2 |  deg_s HH_n^(1) |  dif")
    print("----+-------------+-----------------+-----")

    all_zero = True
    for n in range(1, N + 1):
        L = HH_shift(n, a)
        R = df(n - 1) ** 2 * HH(n, a)
        d = Rs(L - R)
        if d != 0:
            all_zero = False
        print("%3d | %11s | %15s | %s" % (n, df(n - 1) ** 2, L.degree(), d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
