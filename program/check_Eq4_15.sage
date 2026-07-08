# check_Eq4_15.sage -- verification of identity (4.15) (eq:exdf) of d02beta.tex.
#
# For the double factorials a_m = (2m+1)!! the dilated Hankel determinant is
#   HH_n = det( a_{2i+j} )_{0<=i,j<=n-1}
#        = 4^C(n,2) * prod_{i=0}^{n-1} i! (4i+1)!!.
#
# Run:  sage check_Eq4_15.sage N
# prints, for n = 1..N, LHS(n), RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

def df(k):
    # double factorial k!! for odd k >= 1 (empty product = 1)
    return prod(range(k, 0, -2))

def a(m):
    # a_m = (2m+1)!!
    return df(2*m + 1)

def LHS(n):
    # dilated Hankel determinant det(a_{2i+j}), 0<=i,j<=n-1
    return matrix(ZZ, n, n, lambda i, j: a(2*i + j)).det()

def RHS(n):
    return 4^binomial(n, 2) * prod(factorial(i) * df(4*i + 1) for i in range(n))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    print("Identity (4.15) (eq:exdf) of d02beta.tex:  a_m = (2m+1)!!")
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
