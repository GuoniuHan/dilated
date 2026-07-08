# check_Prop29_4.sage -- verification of Proposition 29.4 (prop:hankel-bessel)
# of d57hankel.tex.
#
# Classical Hankel determinant of the Bessel-cosine moments
#   mu_k = (-1)^k (1/2)_k / (s+1)_k       (cosb_s = sum mu_k x^{2k}/(2k)!).
# Proposition 29.4 asserts, for all n>=1 and as an identity in QQ(s),
#   H_n(mu) = det( mu_{i+j} )_{0<=i,j<n}
#           = prod_{i=0}^{n-1} h^S_i
#           = prod_{i=0}^{n-1} (2i)! (s+1/2)_i / ( 4^i (s+i)_i (s+1)_{2i} ).
# First values: H_1 = 1,  H_2 = (2s+1)/(4(s+1)^2(s+2)).
#
# s is kept SYMBOLIC (indeterminate over QQ); H_n is the direct determinant of
# the exact rational moments, checked against the norm-product form and the
# first two tabulated values.
#
# Run:  sage check_Prop29_4.sage N

import sys

Rs = PolynomialRing(QQ, 's'); Fs = Rs.fraction_field()
s = Fs(Rs.gen())

def poch(x, k): return prod(x + j for j in range(k))
H = QQ(1) / 2

def mu(k):
    return (-1) ** k * poch(H, k) / poch(s + 1, k)

def LHS(n):
    return matrix(Fs, n, n, lambda i, j: mu(i + j)).det()

def RHS(n):
    return prod(factorial(2 * i) * poch(s + H, i) / (QQ(4) ** i * poch(s + i, i) * poch(s + 1, 2 * i))
                for i in range(n))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("Proposition 29.4 (prop:hankel-bessel) of d57hankel.tex   [s symbolic]")
    print("H_n(mu) = prod_i (2i)!(s+1/2)_i/(4^i (s+i)_i (s+1)_{2i})   in QQ(s)")
    print("")
    print("  n |  dif = LHS(n) - RHS(n)")
    print("----+-----------------------")

    all_zero = True
    for n in range(1, N + 1):
        d = Fs(LHS(n) - RHS(n))
        if d != 0:
            all_zero = False
        print("%3d | %s" % (n, d))

    # first-value cross-check
    d_h2 = Fs(LHS(2) - (2 * s + 1) / (4 * (s + 1) ** 2 * (s + 2)))
    if d_h2 != 0:
        all_zero = False
    print("")
    print("First-value check  H_2 - (2s+1)/(4(s+1)^2(s+2)) =", d_h2)

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
