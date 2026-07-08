# check_Coro4_5.sage -- verification of Corollary 4.5 (cor:laguerre) of d02beta.tex.
#
# For the Laguerre moments a_m = Gamma(m+alpha+1) the dilated Hankel
# determinant is
#   det( Gamma(2i+j+alpha+1) )_{0<=i,j<=n-1}
#     = Gamma(alpha+1)^n * 2^C(n,2) * prod_{k=1}^{n-1} k!
#         * prod_{m=1}^{n-1} [ (alpha+2m-1)(alpha+2m) ]^{n-m}.
# At a non-negative integer alpha=r this reduces to Corollary 4.4.
#
# alpha is kept exact (symbolic), so the check is a genuine identity, not a
# floating-point approximation.  A non-integer default illustrates the
# continuation.
#
# Run:  sage check_Coro4_5.sage N
# prints, for n = 1..N, LHS(n), RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

# ---------- parameter (any value, integer or not; default 1/2) ----------
alpha = QQ(1)/2
# ------------------------------------------------------------------------

def a(m):
    # a_m = Gamma(m+alpha+1)
    return gamma(m + alpha + 1)

def LHS(n):
    # dilated Hankel determinant det(a_{2i+j}), 0<=i,j<=n-1
    return matrix(SR, n, n, lambda i, j: a(2*i + j)).det()

def RHS(n):
    return (gamma(alpha + 1)^n
            * 2^binomial(n, 2)
            * prod(factorial(k) for k in range(1, n))
            * prod(((alpha + 2*m - 1)*(alpha + 2*m))^(n - m) for m in range(1, n)))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    print("Corollary 4.5 (cor:laguerre) of d02beta.tex")
    print("parameter: alpha = %s" % alpha)
    print("")
    print("  n |          LHS(n)          |          RHS(n)          |  dif")
    print("----+--------------------------+--------------------------+-----")

    all_zero = True
    for n in range(1, N + 1):
        L = SR(LHS(n)).simplify_full()
        R = SR(RHS(n)).simplify_full()
        d = (L - R).simplify_full()
        if d != 0:
            all_zero = False
        print("%3d | %24s | %24s | %s" % (n, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
