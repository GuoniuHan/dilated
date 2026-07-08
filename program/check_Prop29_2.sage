# check_Prop29_2.sage -- verification of Proposition 29.2 (prop:hankel-sineodd)
# of d57hankel.tex.
#
# Classical Hankel determinant of mu_k = (2k+1) b_k, where b_k are the scaled
# Taylor coefficients of the even part x/sin x of the reciprocal sine
# (x/sin x = sum b_k x^{2k}/(2k)!), so
#   mu_0,mu_1,mu_2,mu_3,... = 1, 1, 7/3, 31/3, ...
# Proposition 29.2 asserts, for all n>=1,
#   H_n(mu) = det( (2(i+j)+1) b_{i+j} )_{0<=i,j<n}
#           = prod_{k=1}^{n-1} ( 4k^6/((2k-1)(2k+1)) )^{n-k}
#           = prod_{k=0}^{n-1} 16^k (k!)^8 / ((2k)! (2k+1)!).
# First values: H_1,H_2,H_3,H_4 = 1, 4/3, 4096/135, 50331648/875.
#
# b_k are generated exactly from the x/sin x series; H_n is the direct
# determinant (exact rational arithmetic), checked against both product forms
# and the tabulated first values.
#
# Run:  sage check_Prop29_2.sage N

import sys

def build_b(K):
    R = PowerSeriesRing(QQ, 'x', default_prec=4 * K + 4)
    x = R.gen()
    g = 1 / (sin(x) / x)                      # x / sin x
    return [factorial(2 * k) * g[2 * k] for k in range(K + 1)]

VALS = [1, QQ(4) / 3, QQ(4096) / 135, QQ(50331648) / 875]

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    b = build_b(2 * N)
    mu = [(2 * k + 1) * b[k] for k in range(2 * N)]

    print("Proposition 29.2 (prop:hankel-sineodd) of d57hankel.tex")
    print("H_n( (2k+1)b_k ) = prod (4k^6/((2k-1)(2k+1)))^{n-k}")
    print("                 = prod 16^k (k!)^8/((2k)!(2k+1)!)")
    print("")
    print("  n |               H_n               | dif(form1) | dif(form2) | dif(vals)")
    print("----+---------------------------------+------------+------------+----------")

    all_zero = True
    for n in range(1, N + 1):
        H = matrix(QQ, n, n, lambda i, j: mu[i + j]).det()
        f1 = prod((QQ(4 * k ** 6) / ((2 * k - 1) * (2 * k + 1))) ** (n - k) for k in range(1, n))
        f2 = prod(QQ(16) ** k * factorial(k) ** 8 / (factorial(2 * k) * factorial(2 * k + 1))
                  for k in range(n))
        d1, d2 = H - f1, H - f2
        dv = H - VALS[n - 1] if n - 1 < len(VALS) else QQ(0)
        if d1 or d2 or dv:
            all_zero = False
        print("%3d | %31s | %10s | %10s | %s" % (n, H, d1, d2, dv))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
