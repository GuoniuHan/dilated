# check_Prop24_2.sage -- verification of Proposition 24.2 (prop:besseldshift) of d24Bessel.tex.
#
# The double shift of the Bessel (s,t) family f_{s,t} = cosb_s + int_0^x cosb_t.
# Its coefficients a_m = m! [x^m] f_{s,t} are (eq:bessel-moments)
#   a_{2k}   = (-1)^k (1/2)_k / (s+1)_k,
#   a_{2k+1} = (-1)^k (1/2)_k / (t+1)_k        (k >= 0).
# With nb = ceil(n/2), nn = floor(n/2), the proposition asserts
#   HH_n^{(2)} = det(a_{2i+j+2})_{0<=i,j<n}
#              = (-1)^n (2n-1)!! / ( 2^n (s+n)_nb (t+n)_nn ) * HH_n,
# where HH_n = det(a_{2i+j}) is the unshifted determinant of Theorem 23.1.
#
# Both s and t are kept SYMBOLIC: the moments are rational functions in QQ(s,t),
# and the identity is checked over QQ(s,t).
#
# LHS(n) = HH_n^{(2)};  RHS(n) = the double-factorial factor times HH_n.
#
# Run:  sage check_Prop24_2.sage N
# prints, for n = 1..N, dif = LHS(n)-RHS(n), then a resume.

import sys

Rst = PolynomialRing(QQ, 's, t')
s, t = Rst.gens()
F = Rst.fraction_field()

def poch(a, k):
    return prod(a + j for j in range(k))          # Pochhammer (a)_k

def df(k):
    return prod(range(k, 0, -2))                  # double factorial k!!

def a(m):
    # Bessel moments a_m in QQ(s,t)
    k = m // 2
    if m % 2 == 0:
        return F((-1) ** k * poch(QQ(1)/2, k) / poch(s + 1, k))
    return F((-1) ** k * poch(QQ(1)/2, k) / poch(t + 1, k))

def HH(n):
    return matrix(F, n, n, lambda i, j: a(2*i + j)).det()

def HH_dshift(n):
    return matrix(F, n, n, lambda i, j: a(2*i + j + 2)).det()

def factor(n):
    nb = (n + 1) // 2          # ceil(n/2)
    nn = n // 2                # floor(n/2)
    return (-1) ** n * df(2*n - 1) / (2 ** n * poch(s + n, nb) * poch(t + n, nn))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("Proposition 24.2 (prop:besseldshift) of d24Bessel.tex   [s, t kept symbolic]")
    print("HH_n^(2) = (-1)^n (2n-1)!! / (2^n (s+n)_nb (t+n)_nn) * HH_n")
    print("")
    print("  n |  dif = HH_n^(2) - factor * HH_n")
    print("----+--------------------------------")

    all_zero = True
    for n in range(1, N + 1):
        L = HH_dshift(n)
        R = factor(n) * HH(n)
        d = F(L - R)
        if d != 0:
            all_zero = False
        print("%3d | %s" % (n, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
