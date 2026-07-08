# check_Lem11_2.sage -- verification of Lemma 11.2 (lem:cdh) of d06xcos.tex.
#
# The functional T* at s = 0 is classical (continuous dual Hahn).
# T* is the moment functional in the variable y with
#   T*[y^k] = (2k+1) E_{2k}   (= 1, 3, 25, 427, 12465, ...).
# The lemma asserts that its monic orthogonal polynomials rho_m satisfy
#   rho_{m+1} = (y - 8m^2 - 8m - 3) rho_m - (2m)^4 rho_{m-1},
# with orthogonality and norms
#   T*[rho_m rho_{m'}] = delta_{mm'} ((2m)!!)^4.
#
# We build rho_m from the stated recurrence and verify the orthogonality
# relations and norms against T*.  Everything is exact over QQ.
#
# Indexed by pairs; "sage check_Lem11_2.sage N" checks all 0 <= m,m' <= N.
#
# Run:  sage check_Lem11_2.sage N
# prints, for each pair, LHS, RHS and dif = LHS-RHS,
# then a resume stating whether every dif is 0.

import sys

def df(k):
    return prod(range(k, 0, -2))          # double factorial k!!

def moments_Tstar(M):
    # T*[y^k] = (2k+1) E_{2k},  k = 0..M, exact over QQ
    Rx = PowerSeriesRing(QQ, 'x', default_prec=2 * M + 3)
    x = Rx.gen()
    sec = 1 / cos(x)
    return [(2 * k + 1) * factorial(2 * k) * sec[2 * k] for k in range(M + 1)]

Ry = PolynomialRing(QQ, 'y')
y = Ry.gen()

# recurrence coefficients of rho_m:  c*_m = 8m^2+8m+3,  lambda*_m = (2m)^4
cStar = lambda m: 8*m**2 + 8*m + 3
lStar = lambda m: (2*m)**4

def build_rho(N):
    rho = [Ry(1)]
    if N >= 1:
        rho.append(y - cStar(0))
    for k in range(1, N):
        rho.append((y - cStar(k)) * rho[k] - lStar(k) * rho[k - 1])
    return rho

def apply_T(poly, mu):
    return sum(coef * mu[k] for k, coef in enumerate(poly.list()))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    mu = moments_Tstar(2 * N)
    rho = build_rho(N)

    print("Lemma 11.2 (lem:cdh) of d06xcos.tex   [exact over QQ]")
    print("T*[y^k] = (2k+1)E_{2k} = %s ..." % mu[:min(N + 1, 6)])
    print("")
    print("T*[rho_m rho_m']  vs  delta_{mm'} ((2m)!!)^4")
    print("(m,m') |             LHS             |             RHS             |  dif")
    print("-------+-----------------------------+-----------------------------+-----")

    all_zero = True
    for m in range(N + 1):
        for mp in range(N + 1):
            L = apply_T(rho[m] * rho[mp], mu)
            R = df(2 * m) ** 4 if m == mp else 0
            d = L - R
            if d != 0:
                all_zero = False
            print("(%d,%d)  | %27s | %27s | %s" % (m, mp, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= m,m' <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
