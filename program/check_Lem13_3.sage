# check_Lem13_3.sage -- verification of Lemma 13.3 (lem:oddred) of d07xcosshift.tex.
#
# Odd-sided reduction of the shifted determinant.  For f = (1+x)/cos(x)^{s+1},
# b_m := a_{m+1}, so HH_n^{(1)} = det(b_{2i+j}) with even part
#   b_{2k}   = (2k+1) E^{(s+1)}_{2k} = T*[y^k],
#   b_{2k+1} = E^{(s+1)}_{2k+2}      = S1[y^k]  (S1 = Christoffel transform y.S).
# The lemma asserts, with nb = ceil(n/2), nn = floor(n/2),
#   HH_n^{(1)} = (-1)^C(nn+1,2) * (prod_{m=0}^{nn-1} h^{S1}_m)
#                 * det( T*[ Ph_{nn+r} y^l ] )_{0<=r,l<=nb-1},
# where Ph_m are the monic S1-orthogonal polynomials (kernel polynomials of S)
# and h^{S1}_m = S1[Ph_m^2] their norms.
#
# s is kept SYMBOLIC over QQ.  Ph_m are built here as the kernel polynomials
# Ph_a = (P_{a+1} - nu_a P_a)/y, nu_a = P_{a+1}(0)/P_a(0), from the monic
# S-orthogonal polynomials P_i of eq:Sdata; h^{S1}_m is taken as S1[Ph_m^2].
#
# LHS(n) = HH_n^{(1)} = det(a_{2i+j+1}); RHS(n) = the reduction above.
#
# Run:  sage check_Lem13_3.sage N
# prints, for n = 1..N, dif = LHS(n)-RHS(n), then a resume.

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()
Fs = Rs.fraction_field()
Ry = PolynomialRing(Fs, 'y')
y = Ry.gen()

def build_moments(M):
    # E^{(s+1)}_{2k} = (2k)! [x^{2k}] cos(x)^{-(s+1)},  k = 0..M
    Rx = PowerSeriesRing(Rs, 'x', default_prec=2 * M + 4)
    x = Rx.gen()
    Ssec = (-(s + 1) * cos(x).log()).exp()
    return [factorial(2 * k) * Ssec[2 * k] for k in range(M + 2)]

# monic S-orthogonal polynomials P_i (eq:Sdata)
cS = lambda i: (4*i + 1)*s + 8*i**2 + 4*i + 1
lS = lambda i: 2*i*(2*i - 1)*(s + 2*i - 1)*(s + 2*i)

def buildP(N):
    P = [Ry(1)]
    if N >= 1:
        P.append(y - cS(0))
    for k in range(1, N):
        P.append((y - cS(k)) * P[k] - lS(k) * P[k - 1])
    return P

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    E = build_moments(2 * N + 2)
    Smom  = lambda k: Fs(E[k])            # S[y^k]  = E_{2k}
    S1mom = lambda k: Fs(E[k + 1])        # S1[y^k] = E_{2k+2}
    Tsmom = lambda k: Fs((2*k + 1) * E[k])  # T*[y^k] = (2k+1)E_{2k}
    apply_mom = lambda p, mom: sum(c * mom(k) for k, c in enumerate(p.list()))

    P = buildP(2 * N + 2)
    def Ph(a):                            # kernel polynomial (S1-orthogonal)
        nu = P[a + 1].subs(y=0) / P[a].subs(y=0)
        return (P[a + 1] - nu * P[a]).shift(-1)

    def av(m):                            # a_m
        k = m // 2
        return Fs(E[k]) if m % 2 == 0 else Fs((2*k + 1) * E[k])
    def HH_shift(n):
        return matrix(Fs, n, n, lambda i, j: av(2*i + j + 1)).det()

    print("Lemma 13.3 (lem:oddred) of d07xcosshift.tex   [s kept symbolic]")
    print("")
    print("  n |  dif")
    print("----+-----")

    all_zero = True
    for n in range(1, N + 1):
        nb = (n + 1) // 2
        nn = n // 2
        sign = (-1) ** binomial(nn + 1, 2)
        hnorm = prod(apply_mom(Ph(m)**2, S1mom) for m in range(nn))
        Dmat = matrix(Fs, nb, nb,
                      lambda r, l: apply_mom(Ph(nn + r) * y**l, Tsmom))
        R = sign * hnorm * Dmat.det()
        d = HH_shift(n) - R
        if d != 0:
            all_zero = False
        print("%3d | %s" % (n, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
