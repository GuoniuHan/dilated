# check_Lem11_4.sage -- verification of Lemma 11.4 (lem:conn3) of d06xcos.tex.
#
# Connection formula between the two classical families at s = 0.  Let Ph_i be
# the monic secant (S-orthogonal) polynomials of Lemma 7.2 (lem:cf),
#   Ph_{i+1} = (y - (2i)^2 - (2i+1)^2) Ph_i - ((2i-1)(2i))^2 Ph_{i-1},
# and rho_m the monic T*-orthogonal polynomials of Lemma 11.2 (lem:cdh),
#   rho_{m+1} = (y - 8m^2 - 8m - 3) rho_m - (2m)^4 rho_{m-1}.
# Writing Ph_i = sum_m ktil_{i,m} rho_m, the lemma asserts
#   ktil_{i,m}     = (i!/m!)^2 * 4^{i-m} * binomial(1/2, i-m),
#   T*[Ph_i rho_m] = (i!)^2 (m!)^2 * 2^{2i+2m} * binomial(1/2, i-m).
#
# Two independent checks, both exact over QQ:
#   (A) ktil_{i,m} by expanding Ph_i in the rho-basis  vs the formula;
#   (B) T*[Ph_i rho_m] with T*[y^k] = (2k+1)E_{2k}  vs the formula.
#
# Indexed by pairs; "sage check_Lem11_4.sage N" checks all 0 <= i,m <= N
# (ktil = 0 for m > i, so the vanishing is checked too).
#
# Run:  sage check_Lem11_4.sage N
# prints, for each pair, LHS, RHS and dif = LHS-RHS,
# then a resume stating whether every dif is 0.

import sys

def binom_half(d):
    # binomial(1/2, d) = (1/2)(1/2-1)...(1/2-d+1)/d!  for d >= 0;  0 for d < 0
    if d < 0:
        return QQ(0)
    return prod(QQ(1)/2 - j for j in range(d)) / factorial(d)

def moments_Tstar(M):
    Rx = PowerSeriesRing(QQ, 'x', default_prec=2 * M + 3)
    x = Rx.gen()
    sec = 1 / cos(x)
    return [(2 * k + 1) * factorial(2 * k) * sec[2 * k] for k in range(M + 1)]

Ry = PolynomialRing(QQ, 'y')
y = Ry.gen()

# secant polynomials Ph_i (S, s=0)
cS = lambda i: (2*i)**2 + (2*i + 1)**2
lS = lambda i: ((2*i - 1) * (2*i))**2
# T* polynomials rho_m
cStar = lambda m: 8*m**2 + 8*m + 3
lStar = lambda m: (2*m)**4

def build_family(N, cfun, lfun):
    p = [Ry(1)]
    if N >= 1:
        p.append(y - cfun(0))
    for k in range(1, N):
        p.append((y - cfun(k)) * p[k] - lfun(k) * p[k - 1])
    return p

def expand(Pi, basis, i):
    # coefficients of monic Pi in the monic basis:  Pi = sum_m kap[m] basis[m]
    rem = Pi
    kap = [QQ(0)] * (i + 1)
    for m in range(i, -1, -1):
        cm = rem[m]
        kap[m] = cm
        rem = rem - cm * basis[m]
    return kap

def apply_T(poly, mu):
    return sum(coef * mu[k] for k, coef in enumerate(poly.list()))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    Ph = build_family(N, cS, lS)
    rho = build_family(N, cStar, lStar)
    mu = moments_Tstar(2 * N)
    kexp = [expand(Ph[i], rho, i) for i in range(N + 1)]

    def ktil_formula(i, m):
        return (factorial(i) / factorial(m))**2 * 4**(i - m) * binom_half(i - m)

    print("Lemma 11.4 (lem:conn3) of d06xcos.tex   [exact over QQ]")
    print("")

    all_zero = True

    print("(A) ktil_{i,m}:  expansion of Ph_i in rho-basis  vs  (i!/m!)^2 4^{i-m} C(1/2,i-m)")
    print("(i,m) |             LHS             |             RHS             |  dif")
    print("------+-----------------------------+-----------------------------+-----")
    for i in range(N + 1):
        for m in range(N + 1):
            L = kexp[i][m] if m <= i else QQ(0)
            R = ktil_formula(i, m)
            d = L - R
            if d != 0:
                all_zero = False
            print("(%d,%d) | %27s | %27s | %s" % (i, m, L, R, d))

    print("")
    print("(B) T*[Ph_i rho_m]  vs  (i!)^2 (m!)^2 2^{2i+2m} C(1/2,i-m)")
    print("(i,m) |             LHS             |             RHS             |  dif")
    print("------+-----------------------------+-----------------------------+-----")
    for i in range(N + 1):
        for m in range(N + 1):
            L = apply_T(Ph[i] * rho[m], mu)
            R = factorial(i)**2 * factorial(m)**2 * 2**(2*i + 2*m) * binom_half(i - m)
            d = L - R
            if d != 0:
                all_zero = False
            print("(%d,%d) | %27s | %27s | %s" % (i, m, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= i,m <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
