# check_Lem13_8.sage -- verification of Lemma 13.8 (lem:samekappa) of d07xcosshift.tex.
#
# Same connection coefficients for the half-step T* -> S1 (all at s = 0).
# rho_i are the monic T*-orthogonal polynomials (CDH(1/2,1/2,1/2)),
#   rho_{i+1} = (y - 8i^2 - 8i - 3) rho_i - (2i)^4 rho_{i-1},
# and Ph_m the monic S1-orthogonal kernel polynomials (CDH(1/2,1/2,1)) at s=0,
#   Ph_{m+1} = (y - 8m^2 - 12m - 5) Ph_m - ((2m)(2m+1))^2 Ph_{m-1}.
# Writing rho_i = sum_m kappa^b_{i,m} Ph_m, the lemma asserts
#   kappa^b_{i,m} = ktil_{i,m} := (i!/m!)^2 4^{i-m} binomial(1/2, i-m).
#
# Everything is exact over QQ (s = 0).  kappa^b_{i,m} is obtained by expanding
# rho_i in the Ph-basis and compared with the closed form ktil_{i,m}.
#
# Indexed by pairs; "sage check_Lem13_8.sage N" checks all 0 <= i,m <= N
# (ktil = 0 for m > i, so the vanishing is checked too).
#
# Run:  sage check_Lem13_8.sage N
# prints, for each pair, LHS, RHS and dif = LHS-RHS, then a resume.

import sys

Ry = PolynomialRing(QQ, 'y')
y = Ry.gen()

def binom_half(d):
    if d < 0:
        return QQ(0)
    return prod(QQ(1)/2 - j for j in range(d)) / factorial(d)

# T* polynomials rho_i (s=0)
cT = lambda i: 8*i**2 + 8*i + 3
lT = lambda i: (2*i)**4
# S1 kernel polynomials Ph_m (s=0)
cS1 = lambda m: 8*m**2 + 12*m + 5
lS1 = lambda m: ((2*m) * (2*m + 1))**2

def build_family(N, cfun, lfun):
    p = [Ry(1)]
    if N >= 1:
        p.append(y - cfun(0))
    for k in range(1, N):
        p.append((y - cfun(k)) * p[k] - lfun(k) * p[k - 1])
    return p

def expand(Pi, basis, i):
    rem = Pi
    kap = [QQ(0)] * (i + 1)
    for m in range(i, -1, -1):
        cm = rem[m]
        kap[m] = cm
        rem = rem - cm * basis[m]
    return kap

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    rho = build_family(N, cT, lT)
    Ph = build_family(N, cS1, lS1)
    kexp = [expand(rho[i], Ph, i) for i in range(N + 1)]

    def ktil(i, m):
        return QQ((factorial(i) / factorial(m))**2) * QQ(4)**(i - m) * binom_half(i - m)

    print("Lemma 13.8 (lem:samekappa) of d07xcosshift.tex   [exact over QQ, s=0]")
    print("kappa^b_{i,m} (rho_i in Ph-basis)  vs  (i!/m!)^2 4^{i-m} C(1/2,i-m)")
    print("")
    print("(i,m) |             LHS             |             RHS             |  dif")
    print("------+-----------------------------+-----------------------------+-----")

    all_zero = True
    for i in range(N + 1):
        for m in range(N + 1):
            L = kexp[i][m] if m <= i else QQ(0)
            R = ktil(i, m)
            d = L - R
            if d != 0:
                all_zero = False
            print("(%d,%d) | %27s | %27s | %s" % (i, m, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= i,m <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
