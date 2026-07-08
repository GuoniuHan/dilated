# check_Coro25_6.sage -- verification of Corollary 25.6 (cor:xbe-sinc) of d25xBessel.tex.
#
# The squared-sinc member nu = 1/2:  f = (1+x)(sin x / x)^2.  Its even moments are
#   a_{2k} = (-4)^k (1/2)_k (1)_k / ((3/2)_k (2)_k),   a_{2k+1} = (2k+1) a_{2k}
# (Lemma 25.1 at nu=1/2).  The corollary asserts, for n = 2N,
#   HH_{2N} = (-1)^N 2^{8N^2-5N} (prod_{k=1}^{2N-1} k!)^2 (N-1)!/(3N-1)!
#             * prod_{i=0}^{2N-1} (2i)!/(2N+2i-1)!,
# with the explicit values -2^2/3, 2^16/(3^4 5^3 7^2), -2^47/(3^9 5^5 7^5 11^3 13^2),
# 2^88/(3^14 5^7 7^7 11^5 13^4 17^3 19^2) at n = 2,4,6,8.
#
# All exact over QQ (nu = 1/2 specialised).  Argument = max half-order N.
#
# Run:  sage check_Coro25_6.sage N    (checks half-orders 1..N)

import sys

def poch(base, k):
    return prod(base + j for j in range(k))

def a2k_half(k):
    return (QQ((-4) ** k) * poch(QQ(1) / 2, k) * poch(QQ(1), k)
            / (poch(QQ(3) / 2, k) * poch(QQ(2), k)))

def a(m):
    k = m // 2
    return a2k_half(k) if m % 2 == 0 else (2 * k + 1) * a2k_half(k)

def RHS(N):
    n = 2 * N
    return (QQ((-1) ** N) * QQ(2) ** (8 * N ** 2 - 5 * N)
            * prod(factorial(k) for k in range(1, n)) ** 2
            * factorial(N - 1) / factorial(3 * N - 1)
            * prod(factorial(2 * i) / factorial(2 * N + 2 * i - 1) for i in range(n)))

# explicit values from the paper, indexed by N (order n=2N), N=1..4
LITERAL = {
    1: QQ(-4) / 3,
    2: QQ(2) ** 16 / (QQ(3) ** 4 * QQ(5) ** 3 * QQ(7) ** 2),
    3: QQ(-1) * QQ(2) ** 47 / (QQ(3) ** 9 * QQ(5) ** 5 * QQ(7) ** 5 * QQ(11) ** 3 * QQ(13) ** 2),
    4: QQ(2) ** 88 / (QQ(3) ** 14 * QQ(5) ** 7 * QQ(7) ** 7 * QQ(11) ** 5 * QQ(13) ** 4
                      * QQ(17) ** 3 * QQ(19) ** 2),
}

if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4

    print("Corollary 25.6 (cor:xbe-sinc) of d25xBessel.tex   [exact over QQ, nu=1/2]")
    print("")
    print("  N | order |            HH_{2N}            | dif(vs formula) | dif(vs literal)")
    print("----+-------+------------------------------+-----------------+----------------")

    all_zero = True
    for N in range(1, Nmax + 1):
        n = 2 * N
        L = matrix(QQ, n, n, lambda i, j: a(2 * i + j)).det()
        d1 = L - RHS(N)
        d2 = (L - LITERAL[N]) if N in LITERAL else "N/A"
        if d1 != 0 or (N in LITERAL and L != LITERAL[N]):
            all_zero = False
        print("%3d | %5d | %28s | %15s | %s" % (N, n, L, d1, d2))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for half-orders 1..%d. Formula checked." % Nmax)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
