# check_Eq4_17.sage -- verification of identity (4.17) (eq:exodd) of d02beta.tex.
#
# For a_m = binomial(2m+1,m) the dilated Hankel determinant is
#   HH_n = det( a_{2i+j} )_{0<=i,j<=n-1}
#        = 2^C(n,2) * prod_{i=0}^{n-1} i! (4i+1)! (n+i)!
#                                     / [ (2i+1)! (2i)! (2i+n)! ].
# (The ordinary Hankel determinant, x_i = i, equals 1.)
#
# Run:  sage check_Eq4_17.sage N
# prints, for n = 1..N, LHS(n), RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

def a(m):
    # a_m = binomial(2m+1,m)
    return binomial(2*m + 1, m)

def LHS(n):
    # dilated Hankel determinant det(a_{2i+j}), 0<=i,j<=n-1
    return matrix(ZZ, n, n, lambda i, j: a(2*i + j)).det()

def RHS(n):
    return 2^binomial(n, 2) * prod(
        factorial(i) * factorial(4*i + 1) * factorial(n + i)
        / (factorial(2*i + 1) * factorial(2*i) * factorial(2*i + n))
        for i in range(n))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    print("Identity (4.17) (eq:exodd) of d02beta.tex:  a_m = binomial(2m+1,m)")
    print("")
    print("  n |            LHS(n)            |            RHS(n)            |  dif")
    print("----+-----------------------------+-----------------------------+-----")

    all_zero = True
    for n in range(1, N + 1):
        L = LHS(n)
        R = RHS(n)
        d = L - R
        if d != 0:
            all_zero = False
        print("%3d | %27s | %27s | %s" % (n, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
