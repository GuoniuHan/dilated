# check_Lem15_3.sage -- verification of Lemma 15.3 (lem:cfshift) of d10runkone.tex.
#
# The shifted-secant family (classical).  In the even variable v, the shifted-secant
# functional is  S1[v^k] := E_{2k+2} = (2k+2)! [x^{2k+2}] sec(x)  (secant numbers
# 1,5,61,1385,...).  The lemma asserts its monic orthogonal polynomials P^(1)_l obey
#   P^(1)_{l+1} = (v - (2l+1)^2 - (2l+2)^2) P^(1)_l - ((2l)(2l+1))^2 P^(1)_{l-1},
# with orthogonality and norms
#   S1[P^(1)_l P^(1)_{l'}] = delta_{ll'} ((2l+1)!)^2.
#
# We build P^(1)_l from the stated recurrence and verify orthogonality/norms
# against S1.  Everything is exact over QQ.  Indexed by pairs (l,l').
#
# Run:  sage check_Lem15_3.sage N   checks all 0 <= l,l' <= N.

import sys

def moments_S1(M):
    # S1[v^k] = E_{2k+2} = (2k+2)! [x^{2k+2}] sec x, k = 0..M, exact over QQ
    Rx = PowerSeriesRing(QQ, 'x', default_prec=2 * M + 5)
    x = Rx.gen()
    sec = 1 / cos(x)
    return [factorial(2 * k + 2) * sec[2 * k + 2] for k in range(M + 1)]

Rv = PolynomialRing(QQ, 'v')
v = Rv.gen()

cP = lambda l: (2 * l + 1) ** 2 + (2 * l + 2) ** 2
lP = lambda l: ((2 * l) * (2 * l + 1)) ** 2

def build_P1(N):
    P = [Rv(1)]
    if N >= 1:
        P.append(v - cP(0))
    for l in range(1, N):
        P.append((v - cP(l)) * P[l] - lP(l) * P[l - 1])
    return P

def apply_S1(poly, mu):
    return sum(coef * mu[k] for k, coef in enumerate(poly.list()))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    mu = moments_S1(2 * N)
    P = build_P1(N)

    print("Lemma 15.3 (lem:cfshift) of d10runkone.tex   [exact over QQ]")
    print("S1[v^k] = E_{2k+2} = %s ..." % mu[:min(N + 1, 6)])
    print("")
    print("S1[P^(1)_l P^(1)_l']  vs  delta_{ll'} ((2l+1)!)^2")
    print("(l,l') |             LHS             |             RHS             |  dif")
    print("-------+----------------------------+----------------------------+-----")

    all_zero = True
    for l in range(N + 1):
        for lp in range(N + 1):
            L = apply_S1(P[l] * P[lp], mu)
            Rn = factorial(2 * l + 1) ** 2 if l == lp else 0
            d = L - Rn
            if d != 0:
                all_zero = False
            print("(%d,%d)  | %26s | %26s | %s" % (l, lp, L, Rn, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= l,l' <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
