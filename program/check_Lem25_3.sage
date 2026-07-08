# check_Lem25_3.sage -- verification of Lemma 25.3 (lem:xbe-coeffdet) of d25xBessel.tex.
#
# The coefficient determinant of the Vandermonde reduction (Prop 25.2) for the
# multiplicative Bessel family f_nu = (1+x) cosb_nu^2 at even orders n = 2N.  With
#   Q_{2l}(x)   = (-4)^l (x+1/2)_l (x+nu+1/2)_l (x+nu+1+l)_{N-1-l} (x+2nu+1+l)_{N-1-l},
#   Q_{2m+1}(x) = (2x+2m+1) Q_{2m}(x),   M_N = (m_{kj}),  Q_j(x)=sum_k m_{kj} x^k,
# the lemma asserts
#   det M_N = 2^{N(2N-1)} prod_{k=0}^{N-1} (1/2)_k (nu+1/2)_k^2 (2nu+1/2)_k.
#
# Check: build M_N from the polynomial coefficients and compare det M_N to the closed
# form, as a rational function of the INDETERMINATE nu (dif = 0 in QQ(nu)).
#
# Run:  sage check_Lem25_3.sage N    (checks half-orders 1..N)

import sys

R = PolynomialRing(QQ, 'nu')
nu = R.gen()
K = R.fraction_field()
Kx = PolynomialRing(K, 'x')
x = Kx.gen()

def poch(base, k):
    return prod(base + j for j in range(k))

def build_Q(N):
    Q = [None] * (2 * N)
    for l in range(N):
        Q[2 * l] = Kx(K((-4) ** l) * poch(x + K(1) / 2, l) * poch(x + nu + K(1) / 2, l)
                      * poch(x + nu + 1 + l, N - 1 - l) * poch(x + 2 * nu + 1 + l, N - 1 - l))
    for m in range(N):
        Q[2 * m + 1] = (2 * x + 2 * m + 1) * Q[2 * m]
    return Q

def detM_closed(N):
    return (K(2) ** (N * (2 * N - 1))
            * prod(poch(K(1) / 2, k) * poch(nu + K(1) / 2, k) ** 2 * poch(2 * nu + K(1) / 2, k)
                   for k in range(N)))

if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4

    print("Lemma 25.3 (lem:xbe-coeffdet) of d25xBessel.tex   [exact, nu over QQ(nu)]")
    print("")
    print("  N | order | dif = det M_N - closed form  (should be 0)")
    print("----+-------+--------------------------------------------")

    all_zero = True
    for N in range(1, Nmax + 1):
        n = 2 * N
        Q = build_Q(N)
        MN = matrix(K, n, n, lambda k, j: Q[j][k])
        d = MN.det() - detM_closed(N)
        if d != 0:
            all_zero = False
        print("%3d | %5d | %s" % (N, n, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for half-orders 1..%d. Formula checked." % Nmax)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
