# check_Lem14_3.sage -- verification of Lemma 14.3 (lem:ds-structred) of d08xcosds.tex.
#
# Structure reduction of the residual determinant R_n (double shift of
# f = (1+x)/cos(x)^{s+1}).  With Phat_m the monic S1-orthogonal polynomials
# (data in Lemma 14.2) and the structure form
#   T_i := y Phat_i' - i Phat_i = sum_{a=0}^{i-1} Delta_{i,a} Phat_a,
# Lemma 14.3 asserts, with nb=ceil(n/2), nn=floor(n/2):
#   R_n = 2^{nn} (prod_{a=0}^{nn-1} h^{S1}_a) det( Delta_{nb+r,a} )_{0<=r,a<=nn-1},
# where R_n = det( T1*[Phat_{nb+r} y^m] )_{0<=r,m<=nn-1} is the residual
# determinant of Lemma 14.2.
#
# s is kept SYMBOLIC; each identity is a polynomial identity in s.  Delta_{i,a}
# is computed here directly from its definition (expanding T_i in the Phat
# basis), so this checker is independent of the closed form of Lemma 14.4.
#
# Run:  sage check_Lem14_3.sage N

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

def T1f(p, a): return sum(c * (2 * k + 3) * a[2 * k + 2] for k, c in enumerate(p.list()))

def Delta(i, a_idx, P):
    # coefficient of Phat_{a_idx} in T_i = y Phat_i' - i Phat_i (triangular solve)
    q = y * P[i].derivative() - i * P[i]
    coeffs = {}
    while q != 0 and q.degree() >= 0:
        d = q.degree(); lc = q.leading_coefficient()
        coeffs[d] = lc; q = q - lc * P[d]
    return coeffs.get(a_idx, Rs(0))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(6 * N + 6)
    P = Phat(2 * N + 2)

    print("Lemma 14.3 (lem:ds-structred) of d08xcosds.tex   [s kept symbolic]")
    print("R_n = 2^{nn} (prod h^{S1}_a) det(Delta_{nb+r,a})")
    print("")
    print("  n |  dif = R_n - structure form")
    print("----+-----------------------------")

    all_zero = True
    for n in range(1, N + 1):
        nb = (n + 1) // 2; nn = n // 2
        Rn = matrix(Rs, nn, nn, lambda r, m: T1f(P[nb + r] * y ** m, a)).det()
        Dmat = matrix(Rs, nn, nn, lambda r, ai: Delta(nb + r, ai, P))
        struct = QQ(2) ** nn * prod(hS1(ai) for ai in range(nn)) * Dmat.det()
        d = Rs(Rn - struct)
        if d != 0: all_zero = False
        print("%3d | %s" % (n, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
