# check_Lem9_2.sage -- verification of Lemma 9.2 (lem:mixedconn) of d05eulershift.tex.
#
# Mixed connection at t=1.  Let Q_i be the monic T-orthogonal polynomials for
# v_j = j(j+1) (tangent fraction, t=1), and let P'_m be the monic
# S'-orthogonal polynomials -- the Christoffel transform S' = y*S of the
# secant functional S with u_j = j(j+s) (Lemma lem:christoffel).  Their
# three-term recurrences are
#   Q_{i+1}  = (y - c^T_i)  Q_i  - lam^T_i  Q_{i-1},
#     c^T_i  = 2(2i+1)^2,  lam^T_i = (2i-1)(2i)^2(2i+1);
#   P'_{m+1} = (y - c^{S'}_m) P'_m - lam^{S'}_m P'_{m-1},
#     c^{S'}_m = (2m+1)(2m+1+s) + (2m+2)(2m+2+s),
#     lam^{S'}_m = (2m)(2m+s)(2m+1)(2m+1+s).
# Lemma 9.2 asserts the expansion Q_i = sum_m kappa~_{i,m} P'_m with
#   kappa~_{i,m} = (2i+1)! / (2m+1)! * binom( (s+1)/2 , i-m ).
#
# s is kept SYMBOLIC.  For each i the identity is a polynomial identity in y
# (coefficients in QQ[s]); the checker reports the difference
#   dif(i) = Q_i - sum_{m=0}^{i} kappa~_{i,m} P'_m
# which must be the zero polynomial.
#
# Run:  sage check_Lem9_2.sage N   (checks i = 0..N)

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()
Ry = PolynomialRing(Rs, 'y')
y = Ry.gen()

def Qpolys(imax):
    # monic T-orthogonal polynomials, t=1
    Q = [Ry(1)]
    if imax >= 1:
        Q.append(y - 2 * (2 * 0 + 1) ** 2)          # c^T_0 = 2
    for i in range(1, imax):
        c = 2 * (2 * i + 1) ** 2
        lam = (2 * i - 1) * (2 * i) ** 2 * (2 * i + 1)
        Q.append((y - c) * Q[i] - lam * Q[i - 1])
    return Q

def Ppolys(imax):
    # monic S'-orthogonal polynomials (Christoffel transform of S)
    P = [Ry(1)]
    if imax >= 1:
        c0 = (1) * (1 + s) + (2) * (2 + s)          # c^{S'}_0
        P.append(y - c0)
    for m in range(1, imax):
        c = (2 * m + 1) * (2 * m + 1 + s) + (2 * m + 2) * (2 * m + 2 + s)
        lam = (2 * m) * (2 * m + s) * (2 * m + 1) * (2 * m + 1 + s)
        P.append((y - c) * P[m] - lam * P[m - 1])
    return P

def kappa(i, m):
    d = i - m
    if d < 0:
        return Rs(0)
    binom = prod(((s + 1) / 2 - r) for r in range(d)) / factorial(d)   # binom((s+1)/2, d)
    return factorial(2 * i + 1) / factorial(2 * m + 1) * binom

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    Q = Qpolys(N + 1)
    P = Ppolys(N + 1)

    print("Lemma 9.2 (lem:mixedconn) of d05eulershift.tex   [s kept symbolic, t=1]")
    print("Q_i = sum_m kappa~_{i,m} P'_m,  kappa~_{i,m} = (2i+1)!/(2m+1)! C((s+1)/2, i-m).")
    print("")
    print("  i | deg Q_i |          dif(i) = Q_i - sum kappa~ P'_m")
    print("----+---------+-----------------------------------------")

    all_zero = True
    for i in range(N + 1):
        rhs = sum(kappa(i, m) * P[m] for m in range(i + 1))
        d = Q[i] - rhs
        if d != 0:
            all_zero = False
        print("%3d | %7d | %s" % (i, Q[i].degree(), d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for i = 0..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
