# check_Prop4_2.sage -- verification of Proposition 4.2 (prop:degenerate) of d02beta.tex.
#
# For the degenerate (beta absent) family a_m = a0 * rho^m * (alpha)_m and ANY
# strictly increasing integers 0 <= x_0 < x_1 < ... < x_{n-1},
#   det( a_{x_i+j} )_{0<=i,j<=n-1}
#     = rho^C(n,2) * prod_{0<=i<j<=n-1} (x_j - x_i) * prod_{i=0}^{n-1} a_{x_i}.
#
# The node sequence is supplied by xnode(i); change it to any strictly
# increasing integer sequence.  The default x_i = i^2 + i is non-arithmetic,
# so this really tests the arbitrary-node statement (not just x_i = r i).
#
# Run:  sage check_Prop4_2.sage N
# prints, for n = 1..N, LHS(n), RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

# ---------- parameters (exact rationals) ----------
a0    = QQ(1)
rho   = QQ(3)
alpha = QQ(1)/2

def xnode(i):
    # any strictly increasing integer sequence of the row indices x_i
    return i^2 + i
# --------------------------------------------------

def poch(x, k):
    return prod(x + j for j in range(k))

def a(m):
    # a_m = a0 * rho^m * (alpha)_m
    return a0 * rho^m * poch(alpha, m)

def LHS(n):
    x = [xnode(i) for i in range(n)]
    return matrix(QQ, n, n, lambda i, j: a(x[i] + j)).det()

def RHS(n):
    x = [xnode(i) for i in range(n)]
    vdm = prod(x[j] - x[i] for i in range(n) for j in range(i + 1, n))
    return rho^binomial(n, 2) * vdm * prod(a(x[i]) for i in range(n))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 8

    # sanity: nodes must be strictly increasing non-negative integers
    xs = [xnode(i) for i in range(N)]
    assert all(v in ZZ and v >= 0 for v in xs) and all(xs[i] < xs[i+1] for i in range(N-1)), \
        "xnode must give strictly increasing non-negative integers"

    print("Proposition 4.2 (prop:degenerate) of d02beta.tex")
    print("parameters: a0 = %s, rho = %s, alpha = %s" % (a0, rho, alpha))
    print("nodes x_i = %s ..." % (", ".join(str(xnode(i)) for i in range(min(N, 6)))))
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
