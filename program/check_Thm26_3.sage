# check_Thm26_3.sage -- verification of Theorem 26.3 (thm:conj54double) of d50chapoton.tex.
#
# From the Chapoton--Han matrix M_N to a dilated Hankel determinant.  With
# A_k = rho(x^k), B_k = rho((1+x)x^k) (rho the index-lowering evaluation, see
# check_Thm26_1.sage), let a be the interleaved sequence a_{2k}=A_k, a_{2k+1}=B_k.
# With nb = ceil(N/2), nn = floor(N/2), tau_N = C(nb,2)+C(nn,2), the theorem asserts
#   det M_N = 2^{tau_N} HH_N(a),   HH_N(a) = det(a_{2i+j})_{0<=i,j<N},
# i.e. the Chapoton--Han determinant is, up to a power of 2, the dilated Hankel
# determinant of (1,1,1,2,2,8,10,64,104,...).
#
# Everything is exact over QQ.  LHS = det M_N (via rho); RHS = 2^{tau_N} det(a_{2i+j}).
#
# Run:  sage check_Thm26_3.sage N
# prints, for N = 1..N, dif = det M_N - 2^{tau_N} HH_N(a), then a resume.

import sys

Rx = PolynomialRing(QQ, 'x')
x = Rx.gen()

def index(P):
    return P.degree() + P.valuation()

def N0(P):
    d = index(P)
    num = (x**d + 1) * P(1) - 2 * P
    q, r = num.quo_rem((x - 1)**2)
    assert r == 0, "N0 division not exact"
    return q

def rho(P):
    if P == 0:
        return QQ(0)
    d = index(P)
    for _ in range(d // 2):
        P = N0(P)
    return P[0]

A = lambda k: rho(x**k)
B = lambda k: rho((1 + x) * x**k)

def a_seq(m):
    k = m // 2
    return A(k) if m % 2 == 0 else B(k)

def M(N):
    return matrix(QQ, N, N, lambda i, j: rho(x**i * (1 + x)**j) * 2**(-(j // 2)))

if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("Theorem 26.3 (thm:conj54double) of d50chapoton.tex   [exact over QQ]")
    print("det M_N = 2^{tau_N} HH_N(a),  a_{2k}=A_k, a_{2k+1}=B_k,  tau_N=C(nb,2)+C(nn,2)")
    print("")
    print("  N | tau_N |     det M_N     |    HH_N(a)    |  dif")
    print("----+-------+-----------------+---------------+-----")

    all_zero = True
    for N in range(1, Nmax + 1):
        nb = (N + 1) // 2
        nn = N // 2
        tau = binomial(nb, 2) + binomial(nn, 2)
        L = M(N).det()
        HHa = matrix(QQ, N, N, lambda i, j: a_seq(2*i + j)).det()
        R = 2**tau * HHa
        d = L - R
        if d != 0:
            all_zero = False
        print("%3d | %5d | %15s | %13s | %s" % (N, tau, L, HHa, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for N = 1..%d. Formula checked." % Nmax)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
