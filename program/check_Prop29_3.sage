# check_Prop29_3.sage -- verification of Proposition 29.3 (prop:hankel-sineeven)
# of d57hankel.tex.
#
# Classical Hankel determinant of the reciprocal-sine coefficients themselves,
#   b_k = (2k)! [x^{2k}](x/sin x),
# so b_0,b_1,b_2,b_3,... = 1, 1/3, 7/15, 31/21, ...
# Proposition 29.3 asserts, for all n>=1,
#   H_n(b) = det( b_{i+j} )_{0<=i,j<n}
#          = prod_{i=0}^{n-1} prod_{j=1}^{2i} j^4/((2j-1)(2j+1))
#          = prod_{l=0}^{n-1} 16^l ((2l)!)^6 / ((4l)! (4l+1)!).
# First values: H_1,H_2,H_3,H_4 = 1, 16/45, 65536/55125, 4294967296/18883865.
#
# b_k are generated exactly from the x/sin x series; H_n is the direct
# determinant (exact rational arithmetic), checked against both product forms
# and the tabulated first values.
#
# Run:  sage check_Prop29_3.sage N

import sys

def build_b(K):
    R = PowerSeriesRing(QQ, 'x', default_prec=4 * K + 4)
    x = R.gen()
    g = 1 / (sin(x) / x)                      # x / sin x
    return [factorial(2 * k) * g[2 * k] for k in range(K + 1)]

VALS = [1, QQ(16) / 45, QQ(65536) / 55125, QQ(4294967296) / 18883865]

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    b = build_b(2 * N)

    print("Proposition 29.3 (prop:hankel-sineeven) of d57hankel.tex")
    print("H_n(b) = prod_i prod_{j<=2i} j^4/((2j-1)(2j+1))")
    print("       = prod_l 16^l ((2l)!)^6/((4l)!(4l+1)!)")
    print("")
    print("  n |               H_n               | dif(form1) | dif(form2) | dif(vals)")
    print("----+---------------------------------+------------+------------+----------")

    all_zero = True
    for n in range(1, N + 1):
        H = matrix(QQ, n, n, lambda i, j: b[i + j]).det()
        f1 = prod(prod(QQ(j ** 4) / ((2 * j - 1) * (2 * j + 1)) for j in range(1, 2 * i + 1))
                  for i in range(n))
        f2 = prod(QQ(16) ** l * factorial(2 * l) ** 6 / (factorial(4 * l) * factorial(4 * l + 1))
                  for l in range(n))
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
