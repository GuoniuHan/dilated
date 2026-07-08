# check_Coro4_7.sage -- verification of Corollary 4.7 (cor:cbinom) of d02beta.tex.
#
# For the central binomial coefficients b_m = binomial(2m,m) the r-step Hankel
# determinant is (eq:cbinomrfold)
#   det( b_{ri+j} )_{0<=i,j<=n-1}
#     = 2^{n-1} * r^C(n,2) * prod_{i=0}^{n-1} i! (2ri)! (n-1+i)!
#                                            / [ (2i)! (ri)! (ri+n-1)! ].
# r = 1 gives the classical det(b_{i+j}) = 2^{n-1}; r = 2 the dilated Hankel
# determinant HH_n(b) = 1,8,224,22528,8200192,... (eq:cbinomdouble).
#
# Run:  sage check_Coro4_7.sage N
# prints, for n = 1..N, LHS(n), RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

# ---------- parameter (step, positive integer; default 2 = dilated) ----------
r = 2
# ----------------------------------------------------------------------------

def b(m):
    # central binomial coefficient b_m
    return binomial(2*m, m)

def LHS(n):
    # r-step Hankel determinant det(b_{ri+j}), 0<=i,j<=n-1
    return matrix(ZZ, n, n, lambda i, j: b(r*i + j)).det()

def RHS(n):
    return (2^(n - 1) * r^binomial(n, 2)
            * prod(factorial(i) * factorial(2*r*i) * factorial(n - 1 + i)
                   / (factorial(2*i) * factorial(r*i) * factorial(r*i + n - 1))
                   for i in range(n)))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    print("Corollary 4.7 (cor:cbinom) of d02beta.tex")
    print("parameter: r = %s" % r)
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
