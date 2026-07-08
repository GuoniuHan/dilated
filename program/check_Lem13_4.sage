# check_Lem13_4.sage -- verification of Lemma 13.4 (lem:div) of d07xcosshift.tex.
#
# Divisibility of the shifted determinant.  For f = (1+x)/cos(x)^{s+1} with
#   a_{2k} = E^{(s+1)}_{2k},  a_{2k+1} = (2k+1) E^{(s+1)}_{2k},
# and HH_n^{(1)} = det(a_{2i+j+1})_{0<=i,j<n}, the lemma asserts
#   prod_{j=1}^{n-1} (s+j)^{n-j}  divides  HH_n^{(1)}.
#
# s is kept SYMBOLIC over QQ.  We check divisibility exactly: the remainder of
# HH_n^{(1)} on division by D_n = prod_{j=1}^{n-1}(s+j)^{n-j} is 0, and the
# quotient HH_n^{(1)}/D_n is a genuine polynomial in s (its printed here).
#
# Run:  sage check_Lem13_4.sage N
# prints, for n = 1..N, the remainder (dif) and the quotient, then a resume.

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()

def build_a(M):
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2)
    x = Rx.gen()
    g = (1 + x) * (-(s + 1) * cos(x).log()).exp()
    return [factorial(m) * g[m] for m in range(M + 1)]

def HH_shift(n, a):
    return matrix(Rs, n, n, lambda i, j: a[2 * i + j + 1]).det()

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)

    print("Lemma 13.4 (lem:div) of d07xcosshift.tex   [s kept symbolic]")
    print("prod_{j=1}^{n-1} (s+j)^{n-j}  divides  HH_n^(1)")
    print("")
    print("  n |               quotient HH_n^(1)/D_n               | remainder (dif)")
    print("----+---------------------------------------------------+----------------")

    all_zero = True
    for n in range(1, N + 1):
        H = HH_shift(n, a)
        D = prod((s + j) ** (n - j) for j in range(1, n))
        q, r = H.quo_rem(Rs(D))
        if r != 0:
            all_zero = False
        print("%3d | %49s | %s" % (n, q, r))

    print("")
    if all_zero:
        print("Resume: all remainders = 0 for n = 1..%d. Divisibility checked." % N)
    else:
        print("Resume: some remainder != 0. Divisibility NOT checked!")
