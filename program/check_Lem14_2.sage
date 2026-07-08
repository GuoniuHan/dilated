# check_Lem14_2.sage -- verification of Lemma 14.2 (lem:ds-red) of d08xcosds.tex.
#
# One-sided reduction of the double shift against the even functional S1.
# For f = (1+x)/cos(x)^{s+1}, put d_m = a_{m+2}, so HH_n^{(2)} = det(d_{2i+j}).
# The even/odd parts of d are the Christoffel transforms
#   S1[y^k]  = a_{2k+2},
#   T1*[y^k] = a_{2k+3} = (2k+3) a_{2k+2},
# and Phat_m are the monic S1-orthogonal polynomials, three-term recurrence
#   Phat_{m+1} = (y - c^{S1}_m) Phat_m - lam^{S1}_m Phat_{m-1},
#   c^{S1}_m = (4m+3)s + 8m^2 + 12m + 5,
#   lam^{S1}_m = 2m(2m+1)(s+2m)(s+2m+1),   h^{S1}_m = (2m+1)! (s+1)_{2m+1}.
# Lemma 14.2 asserts, with nb=ceil(n/2), nn=floor(n/2):
#   (eq:ds-reduction)  HH_n^{(2)} = (-1)^{C(nb,2)} (prod_{l=0}^{nb-1} h^{S1}_l) R_n,
#                      R_n = det( T1*[Phat_{nb+r} y^m] )_{0<=r,m<=nn-1};
#   (eq:ds-entry)      for i>m,  T1*[Phat_i y^m] = 2 S1[y^{m+1} Phat_i'].
#
# s is kept SYMBOLIC; the reduction is a polynomial identity in s (all s),
# the entry identity is checked over all pairs i>m up to N.
#
# Run:  sage check_Lem14_2.sage N

import sys

Rs = PolynomialRing(QQ, 's'); s = Rs.gen()
Ry = PolynomialRing(Rs, 'y'); y = Ry.gen()

def poch(x, k): return prod(x + j for j in range(k))

def build_a(M):
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2); x = Rx.gen()
    logcos = cos(x).log(); Ssec = (-(s + 1) * logcos).exp()
    E = [factorial(2 * k) * Ssec[2 * k] for k in range(M // 2 + 1)]
    a = []
    for m in range(M + 1):
        a.append(E[m // 2] if m % 2 == 0 else m * E[(m - 1) // 2])
    return a

def cS1(m):   return (4 * m + 3) * s + 8 * m ** 2 + 12 * m + 5
def lamS1(m): return 2 * m * (2 * m + 1) * (s + 2 * m) * (s + 2 * m + 1)
def hS1(m):   return factorial(2 * m + 1) * poch(s + 1, 2 * m + 1)

def Phat(mmax):
    P = [Ry(1)]
    if mmax >= 1: P.append(y - cS1(0))
    for m in range(1, mmax):
        P.append((y - cS1(m)) * P[m] - lamS1(m) * P[m - 1])
    return P

def S1f(p, a): return sum(c * a[2 * k + 2] for k, c in enumerate(p.list()))
def T1f(p, a): return sum(c * (2 * k + 3) * a[2 * k + 2] for k, c in enumerate(p.list()))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(6 * N + 6)
    P = Phat(2 * N + 2)

    print("Lemma 14.2 (lem:ds-red) of d08xcosds.tex   [s kept symbolic]")
    print("")
    print("Part 1  (eq:ds-reduction):  HH_n^{(2)} = (-1)^C(nb,2) prod h^{S1}_l * R_n")
    print("  n |  dif = HH_n^{(2)} - reduction")
    print("----+------------------------------")
    all_zero = True
    for n in range(1, N + 1):
        nb = (n + 1) // 2; nn = n // 2
        HH2 = matrix(Rs, n, n, lambda i, j: a[2 * i + j + 2]).det()
        Rn = matrix(Rs, nn, nn, lambda r, m: T1f(P[nb + r] * y ** m, a)).det()
        red = (-1) ** binomial(nb, 2) * prod(hS1(l) for l in range(nb)) * Rn
        d = Rs(HH2 - red)
        if d != 0: all_zero = False
        print("%3d | %s" % (n, d))

    print("")
    print("Part 2  (eq:ds-entry, i>m):  T1*[Phat_i y^m] - 2 S1[y^{m+1} Phat_i']")
    print(" (i,m) |  dif")
    print("-------+-----")
    for i in range(1, N + 1):
        for m in range(i):
            d = Rs(T1f(P[i] * y ** m, a) - 2 * S1f(y ** (m + 1) * P[i].derivative(), a))
            if d != 0: all_zero = False
            print("(%d,%d) | %s" % (i, m, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n,pairs up to %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
