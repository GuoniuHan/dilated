# check_Thm25_4.sage -- verification of Theorem 25.4 (thm:xbesseleven) of d25xBessel.tex.
#
# Closed form of the dilated Hankel determinant of the multiplicative Bessel family
# f_nu = (1+x) cosb_nu^2 at even orders n = 2N.  With
#   a_{2k} = (-4)^k (1/2)_k (nu+1/2)_k / ((nu+1)_k (2nu+1)_k),  a_{2k+1}=(2k+1)a_{2k},
# the theorem asserts, as an identity in QQ(nu),
#   HH_{2N}(f_nu) = (-1)^N 2^{C(2N,2)} prod_{i=0}^{2N-1} (2i)!
#                   * prod_{k=0}^{N-1} (1/2)_k (nu+1/2)_k^2 (2nu+1/2)_k
#                   * prod_{i=0}^{2N-1} (nu+1/2)_i / ((nu+1)_{N-1+i} (2nu+1)_{N-1+i}).
#
# STRONGEST check: nu kept as an INDETERMINATE; LHS = det(a_{2i+j}) and RHS compared
# as rational functions in QQ(nu) (dif = 0).  Also cross-checks the two stated first
# values HH_2 = -2/(nu+1) and HH_4.
#
# Run:  sage check_Thm25_4.sage N    (checks half-orders 1..N)

import sys

R = PolynomialRing(QQ, 'nu')
nu = R.gen()
K = R.fraction_field()

def poch(base, k):
    return prod(base + j for j in range(k))

def a2k(k):
    return (K((-4) ** k) * poch(K(1) / 2, k) * poch(nu + K(1) / 2, k)
            / (poch(nu + 1, k) * poch(2 * nu + 1, k)))

def a(m):
    k = m // 2
    return a2k(k) if m % 2 == 0 else (2 * k + 1) * a2k(k)

def RHS(N):
    n = 2 * N
    return (K((-1) ** N) * K(2) ** binomial(n, 2) * prod(factorial(2 * i) for i in range(n))
            * prod(poch(K(1) / 2, k) * poch(nu + K(1) / 2, k) ** 2 * poch(2 * nu + K(1) / 2, k)
                   for k in range(N))
            * prod(poch(nu + K(1) / 2, i) / (poch(nu + 1, N - 1 + i) * poch(2 * nu + 1, N - 1 + i))
                   for i in range(n)))

if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4

    print("Theorem 25.4 (thm:xbesseleven) of d25xBessel.tex   [exact, nu over QQ(nu)]")
    print("")
    print("  N | order | dif = LHS - RHS  (should be 0)")
    print("----+-------+-------------------------------")

    all_zero = True
    HH = {}
    for N in range(1, Nmax + 1):
        n = 2 * N
        L = matrix(K, n, n, lambda i, j: a(2 * i + j)).det()
        HH[N] = L
        d = L - RHS(N)
        if d != 0:
            all_zero = False
        print("%3d | %5d | %s" % (N, n, d))

    print("")
    # cross-check the two displayed first values
    v2 = K(-2) / (nu + 1)
    v4 = K(2160) * (nu + K(1) / 4) * (nu + K(1) / 2) * (nu + K(5) / 2) \
        / ((nu + 1) ** 7 * (nu + 2) ** 4 * (nu + 3) ** 2 * (nu + 4))
    print("first-value checks:")
    print("  HH_2 - (-2/(nu+1))      = %s" % (HH[1] - v2 if 1 in HH else "N/A"))
    if 2 in HH:
        print("  HH_4 - (2160(...)/(...)) = %s" % (HH[2] - v4))
        if HH[2] - v4 != 0:
            all_zero = False
    if 1 in HH and HH[1] - v2 != 0:
        all_zero = False

    print("")
    if all_zero:
        print("Resume: all dif = 0 for half-orders 1..%d. Formula checked." % Nmax)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
