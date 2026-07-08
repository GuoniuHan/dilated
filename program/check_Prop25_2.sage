# check_Prop25_2.sage -- verification of Proposition 25.2 (prop:xbe-vdm) of d25xBessel.tex.
#
# Row factorisation and Vandermonde reduction of the multiplicative Bessel family
# f_nu = (1+x) cosb_nu^2 at even orders n = 2N.  With
#   a_{2k}   = (-4)^k (1/2)_k (nu+1/2)_k / ((nu+1)_k (2nu+1)_k),  a_{2k+1}=(2k+1)a_{2k},
#   w_i = a_{2i} / ((nu+1+i)_{N-1} (2nu+1+i)_{N-1}),
#   Q_{2l}(x) = (-4)^l (x+1/2)_l (x+nu+1/2)_l (x+nu+1+l)_{N-1-l} (x+2nu+1+l)_{N-1-l},
#   Q_{2m+1}(x) = (2x+2m+1) Q_{2m}(x),
# the proposition asserts (0 <= i,j <= 2N-1):
#   (A) a_{2i+j} = w_i Q_j(i),                                              [eq:xbe-rowfact]
#   (B) HH_{2N}(f_nu) = det M_N * prod_{k=1}^{2N-1} k! * prod_{i=0}^{2N-1} w_i,  [eq:xbe-vdmred]
#       where M_N = (m_{kj}) is the coefficient matrix, Q_j(x) = sum_k m_{kj} x^k.
#
# Everything is exact over QQ(nu) (nu indeterminate).  Argument = max half-order N.
#
# Run:  sage check_Prop25_2.sage N    (checks half-orders 1..N, i.e. orders 2,4,...,2N)

import sys

R = PolynomialRing(QQ, 'nu')
nu = R.gen()
K = R.fraction_field()
Kx = PolynomialRing(K, 'x')
x = Kx.gen()

def poch(base, k):
    return prod(base + j for j in range(k))

def a2k(k):
    return (K((-4) ** k) * poch(K(1) / 2, k) * poch(nu + K(1) / 2, k)
            / (poch(nu + 1, k) * poch(2 * nu + 1, k)))

def a(m):
    k = m // 2
    return a2k(k) if m % 2 == 0 else (2 * k + 1) * a2k(k)

def build_Q(N):
    Q = [None] * (2 * N)
    for l in range(N):
        Q[2 * l] = Kx(K((-4) ** l) * poch(x + K(1) / 2, l) * poch(x + nu + K(1) / 2, l)
                      * poch(x + nu + 1 + l, N - 1 - l) * poch(x + 2 * nu + 1 + l, N - 1 - l))
    for m in range(N):
        Q[2 * m + 1] = (2 * x + 2 * m + 1) * Q[2 * m]
    return Q

def w(i, N):
    return a2k(i) / (poch(nu + 1 + i, N - 1) * poch(2 * nu + 1 + i, N - 1))

if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4

    print("Proposition 25.2 (prop:xbe-vdm) of d25xBessel.tex   [exact, nu over QQ(nu)]")
    print("")
    print("  N | order | (A) rowfact dif | (B) vdmred dif  (both should be 0)")
    print("----+-------+-----------------+----------------------------------")

    all_zero = True
    for N in range(1, Nmax + 1):
        n = 2 * N
        Q = build_Q(N)
        # (A) a_{2i+j} = w_i Q_j(i)
        dA = K(0)
        for i in range(n):
            wi = w(i, N)
            for j in range(n):
                dA += (a(2 * i + j) - wi * Q[j](i)) ** 2
        # (B) reduction
        HH = matrix(K, n, n, lambda i, j: a(2 * i + j)).det()
        MN = matrix(K, n, n, lambda k, j: Q[j][k])       # m_{kj} = [x^k] Q_j
        rhs = MN.det() * prod(factorial(k) for k in range(1, n)) * prod(w(i, N) for i in range(n))
        dB = HH - rhs
        if dA != 0 or dB != 0:
            all_zero = False
        print("%3d | %5d | %15s | %s" % (N, n, dA, dB))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for half-orders 1..%d. Formula checked." % Nmax)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
