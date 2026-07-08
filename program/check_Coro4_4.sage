# check_Coro4_4.sage -- verification of Corollary 4.4 (cor:factorial) of d02beta.tex.
#
# For the factorial sequence a_m = (m+r)! the dilated Hankel determinant is
#       HH_n((m+r)!) = det( (2i+j+r)! )_{0<=i,j<=n-1}
#                    = 2^C(n,2) * prod_{k=1}^{n-1} k! * prod_{i=0}^{n-1} (2i+r)!.
#
# Run:  sage check_Coro4_4.sage N
# prints, for n = 1..N, LHS(n), RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

# ---------- parameter (non-negative integer) ----------
r = 0
# ------------------------------------------------------

def a(m):
    # a_m = (m+r)!
    return factorial(m + r)

def LHS(n):
    # dilated Hankel determinant det(a_{2i+j}), 0<=i,j<=n-1
    return matrix(ZZ, n, n, lambda i, j: a(2*i + j)).det()

def RHS(n):
    return (2^binomial(n, 2)
            * prod(factorial(k) for k in range(1, n))
            * prod(factorial(2*i + r) for i in range(n)))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    print("Corollary 4.4 (cor:factorial) of d02beta.tex")
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
