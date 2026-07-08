# check_Thm6_1.sage -- verification of Theorem 6.1 (thm:dergauss) of d03invo.tex.
#
# Single shift of the Gaussian family.  For f(x) = e^{cx + x^2/2}, with
#   a_n = sum_j binomial(n,2j) (2j-1)!! c^{n-2j},
# the SHIFTED dilated Hankel determinant (that of f') is
#   HH_n(f') = det( a_{2i+j+1} )_{0<=i,j<=n-1}
#            = c^n * HH_n(f) = c^n * (2c)^C(n,2) * prod_{k=1}^{n-1} k!.
# (For c = 1 this equals the unshifted HH_n: invariance under differentiation.)
#
# c is kept SYMBOLIC, so each LHS(n)=RHS(n) is a polynomial identity in c.
#
# Run:  sage check_Thm6_1.sage N
# prints, for n = 1..N, LHS(n), RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

Rc = PolynomialRing(QQ, 'c')
c = Rc.gen()

def df(k):
    # double factorial k!! for odd k (empty product = 1, so (-1)!! = 1)
    return prod(range(k, 0, -2))

def a(n):
    # a_n = sum_j binomial(n,2j) (2j-1)!! c^{n-2j}
    return sum(binomial(n, 2*j) * df(2*j - 1) * c^(n - 2*j)
               for j in range(n//2 + 1))

def LHS(n):
    # shifted dilated Hankel determinant det(a_{2i+j+1})
    return matrix(Rc, n, n, lambda i, j: a(2*i + j + 1)).det()

def RHS(n):
    return c^n * (2*c)^binomial(n, 2) * prod(factorial(k) for k in range(1, n))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    print("Theorem 6.1 (thm:dergauss) of d03invo.tex   [c kept symbolic]")
    print("")
    print("  n |                 LHS(n)                 |                 RHS(n)                 |  dif")
    print("----+----------------------------------------+----------------------------------------+-----")

    all_zero = True
    for n in range(1, N + 1):
        L = LHS(n)
        R = RHS(n)
        d = L - R
        if d != 0:
            all_zero = False
        print("%3d | %38s | %38s | %s" % (n, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
