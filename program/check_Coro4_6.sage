# check_Coro4_6.sage -- verification of Corollary 4.6 (cor:catalan) of d02beta.tex.
#
# For the Catalan numbers C_m = binomial(2m,m)/(m+1) the r-step Hankel
# determinant is (eq:catalanrfold)
#   det( C_{ri+j} )_{0<=i,j<=n-1}
#     = r^C(n,2) * prod_{i=0}^{n-1} i! (2ri)! (n+i)! / [ (2i)! (ri)! (ri+n)! ].
# r = 1 gives the classical det(C_{i+j}) = 1; r = 2 the dilated Hankel
# determinant HH_n(C) = 1,3,32,1232,172032,... (eq:catalandouble);
# r = 3 the 3-step values 1,9,990,1363230,...
#
# Run:  sage check_Coro4_6.sage N
# prints, for n = 1..N, LHS(n), RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

# ---------- parameter (step, positive integer; default 2 = dilated) ----------
r = 2
# ----------------------------------------------------------------------------

def C(m):
    # Catalan number C_m
    return binomial(2*m, m) / (m + 1)

def LHS(n):
    # r-step Hankel determinant det(C_{ri+j}), 0<=i,j<=n-1
    return matrix(ZZ, n, n, lambda i, j: C(r*i + j)).det()

def RHS(n):
    return (r^binomial(n, 2)
            * prod(factorial(i) * factorial(2*r*i) * factorial(n + i)
                   / (factorial(2*i) * factorial(r*i) * factorial(r*i + n))
                   for i in range(n)))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    print("Corollary 4.6 (cor:catalan) of d02beta.tex")
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
