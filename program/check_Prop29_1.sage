# check_Prop29_1.sage -- verification of Proposition 29.1 (prop:hankel-secodd)
# of d57hankel.tex.
#
# Classical Hankel determinant of mu_k = (2k+1) E_{2k}, where E_{2k} are the
# secant numbers (1/cos x = sum E_{2k} x^{2k}/(2k)!), so
#   mu_0,mu_1,mu_2,... = 1, 3, 25, 427, 12465, ...
# Proposition 29.1 asserts, for all n>=1,
#   H_n(mu) = det( (2(i+j)+1) E_{2(i+j)} )_{0<=i,j<n}
#           = prod_{k=1}^{n-1} (2k)^{4(n-k)}
#           = 2^{4 C(n,2)} prod_{k=1}^{n-1} (k!)^4.
# First values: H_1,H_2,H_3,H_4 = 1, 2^4, 2^16, 2^32 3^4.
#
# The E_{2k} are generated exactly from the secant series; H_n is the direct
# determinant (exact integer arithmetic), checked against both product forms
# and the tabulated first values.
#
# Run:  sage check_Prop29_1.sage N

import sys

def build_E(K):
    R = PowerSeriesRing(QQ, 'x', default_prec=4 * K + 4)
    x = R.gen()
    sec = 1 / cos(x)
    return [factorial(2 * k) * sec[2 * k] for k in range(K + 1)]

VALS = [1, 2 ** 4, 2 ** 16, 2 ** 32 * 3 ** 4]

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    E = build_E(2 * N)
    mu = [(2 * k + 1) * E[k] for k in range(2 * N)]

    print("Proposition 29.1 (prop:hankel-secodd) of d57hankel.tex")
    print("H_n( (2k+1)E_{2k} ) = prod (2k)^{4(n-k)} = 2^{4C(n,2)} prod (k!)^4")
    print("")
    print("  n |            H_n            | dif(form1) | dif(form2) | dif(vals)")
    print("----+---------------------------+------------+------------+----------")

    all_zero = True
    for n in range(1, N + 1):
        H = matrix(ZZ, n, n, lambda i, j: mu[i + j]).det()
        f1 = prod((2 * k) ** (4 * (n - k)) for k in range(1, n))
        f2 = 2 ** (4 * binomial(n, 2)) * prod(factorial(k) ** 4 for k in range(1, n))
        d1, d2 = H - f1, H - f2
        dv = H - VALS[n - 1] if n - 1 < len(VALS) else 0
        if d1 or d2 or dv:
            all_zero = False
        print("%3d | %25s | %10s | %10s | %s" % (n, H, d1, d2, dv))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
