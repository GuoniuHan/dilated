# check_Thm26_1.sage -- verification of Theorem 26.1 (thm:rootlink) of d50chapoton.tex.
#
# The Chapoton--Han root conjecture equals a dilated Hankel determinant.
# The evaluation rho is Chapoton--Han's index-lowering map (graded, NOT a linear
# functional on QQ[x]): with index(P) = deg P + val P and
#   N0(P) = ((x^d+1) P(1) - 2 P(x)) / (x-1)^2,   d = index(P),
# one sets rho(P) = constant term of N0^{floor(d/2)}(P).  The matrix is
#   M_N(i,j) = rho( x^i (1+x)^j ) * 2^{-floor(j/2)},   0 <= i,j <= N-1.
# The theorem asserts
#   det M_N = 2^{-C(N,2)} HH_N((1+x)/cos x)
#           = ((N-1)!!)^2 prod_{k=1}^{N-2} (k!!)^6
#           = prod_{k=1}^{N-1} ((N-k)!)^{eps(k)},   eps(k)=2 (k odd), 4 (k even).
#
# Everything is exact over QQ.  We check the full chain of equalities:
#   det M_N  vs  the two closed products  vs  2^{-C(N,2)} HH_N((1+x)/cos x)
# (the last computed from the secant sequence a_{2k}=E_{2k}, a_{2k+1}=(2k+1)E_{2k}).
#
# Run:  sage check_Thm26_1.sage N
# prints, for N = 1..N, det M_N and the three difs, then a resume.

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
    return P[0]                      # constant term (index <= 1 here)

def df(k):
    return prod(range(k, 0, -2))     # double factorial k!!

def eps(k):
    return 2 if k % 2 == 1 else 4

def M(N):
    return matrix(QQ, N, N, lambda i, j: rho(x**i * (1 + x)**j) * 2**(-(j // 2)))

if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    # secant numbers E_{2j} and the sequence of (1+x)/cos x
    Rt = PowerSeriesRing(QQ, 'tt', default_prec=2 * Nmax + 3)
    tt = Rt.gen()
    sec = 1 / cos(tt)
    E = lambda j: factorial(2 * j) * sec[2 * j]
    aS = lambda m: E(m // 2) if m % 2 == 0 else (2 * (m // 2) + 1) * E(m // 2)

    print("Theorem 26.1 (thm:rootlink) of d50chapoton.tex   [exact over QQ]")
    print("det M_N = 2^{-C(N,2)} HH_N((1+x)/cos x) = ((N-1)!!)^2 prod (k!!)^6 = CH-product")
    print("")
    print("  N |     det M_N     | dif vs ((N-1)!!)^2.. | dif vs CH-prod | dif vs 2^-C HH")
    print("----+-----------------+----------------------+----------------+---------------")

    all_zero = True
    for N in range(1, Nmax + 1):
        L = M(N).det()
        R1 = df(N - 1)**2 * prod(df(k)**6 for k in range(1, N - 1))
        R2 = prod(factorial(N - k)**eps(k) for k in range(1, N))
        HH = matrix(QQ, N, N, lambda i, j: aS(2*i + j)).det()
        R3 = 2**(-binomial(N, 2)) * HH
        d1, d2, d3 = L - R1, L - R2, L - R3
        if d1 != 0 or d2 != 0 or d3 != 0:
            all_zero = False
        print("%3d | %15s | %20s | %14s | %s" % (N, L, d1, d2, d3))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for N = 1..%d. Formula checked." % Nmax)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
