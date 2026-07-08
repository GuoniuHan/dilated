# check_Coro4_3.sage -- verification of Corollary 4.3 (cor:betadouble) of d02beta.tex.
#
# For the Beta family
#       a_m = a0 * rho^m * (alpha)_m / (beta)_m,
# the dilated Hankel determinant is
#       HH_n(a) = det( a_{2i+j} )_{0<=i,j<=n-1}
#               = (2 rho)^C(n,2) * prod_{i=0}^{n-1} i! (beta-alpha)_i a_{2i}/(2i+beta)_{n-1}.
#
# Run:  sage check_Coro4_3.sage N
# prints, for n = 1..N, LHS(n), RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

# ---------- parameters (exact rationals; beta not in {0,-1,-2,...}) ----------
a0    = QQ(1)
rho   = QQ(3)
alpha = QQ(1)/2
beta  = QQ(2)
# ---------------------------------------------------------------------------

def poch(x, k):
    return prod(x + j for j in range(k))

def a(m):
    # a_m = a0 * rho^m * (alpha)_m / (beta)_m
    return a0 * rho^m * poch(alpha, m) / poch(beta, m)

def LHS(n):
    # dilated Hankel determinant det(a_{2i+j}), 0<=i,j<=n-1
    return matrix(QQ, n, n, lambda i, j: a(2*i + j)).det()

def RHS(n):
    return (2*rho)^binomial(n, 2) * prod(
        factorial(i) * poch(beta - alpha, i) * a(2*i) / poch(2*i + beta, n - 1)
        for i in range(n)
    )

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    print("Corollary 4.3 (cor:betadouble) of d02beta.tex")
    print("parameters: a0 = %s, rho = %s, alpha = %s, beta = %s"
          % (a0, rho, alpha, beta))
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
