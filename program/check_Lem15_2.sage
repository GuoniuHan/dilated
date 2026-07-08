# check_Lem15_2.sage -- verification of Lemma 15.2 (lem:rank1) of d10runkone.tex.
#
# Rank-one structure of the perturbation f_s = (sin x+1)/cos^2 x + s sin x.
# With a_n = E_{n+1} (shifted Euler), A = (a_{2i+j})_{0<=i,j<=n-1}, and
#   sigma_{2k}=0, sigma_{2k+1}=(-1)^k  (coefficients of sin x),
# the lemma asserts:
#   (A) sigma_{2i+j} = (-1)^i sigma_j        for all i,j;
#   (B) the perturbation matrix (sigma_{2i+j})_{i,j} = u w^T,  u_i=(-1)^i, w_j=sigma_j;
#   (C) matrix determinant lemma:
#         HH_n(f_s) = det(A + s u w^T) = det(A)(1 + s w^T A^{-1} u),
#       with det(A) = HH_n^E = prod_{k=1}^{n-1}(k!)^2(2k+1)!!.
#
# Checks (A),(B) are exact identities over QQ; (C) compares det(A + s u w^T) with
# det(A)(1 + s w^T A^{-1} u) as polynomials in s over QQ (dif is the zero poly).
#
# Run:  sage check_Lem15_2.sage N
# prints (A) over pairs (i,j), (B),(C) over n=1..N, then a resume.

import sys

R = PolynomialRing(QQ, 's')
s = R.gen()

def df(k):
    return prod(range(k, 0, -2))

def build_a(M):
    Rx = PowerSeriesRing(QQ, 'x', default_prec=M + 2)
    x = Rx.gen()
    g = (sin(x) + 1) / cos(x) ** 2
    return [factorial(m) * g[m] for m in range(M + 1)]

def sigma(m):
    if m % 2 == 0:
        return QQ(0)
    return QQ((-1) ** ((m - 1) // 2))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)

    def HE(n):
        return prod(factorial(k) ** 2 * df(2 * k + 1) for k in range(1, n))

    all_ok = True

    print("Lemma 15.2 (lem:rank1) of d10runkone.tex   [exact over QQ / QQ[s]]")
    print("")
    print("(A) sigma_{2i+j} = (-1)^i sigma_j    for 0 <= i,j <= N")
    print("(i,j) |  LHS  |  RHS  | dif")
    print("------+-------+-------+----")
    for i in range(N + 1):
        for j in range(N + 1):
            L = sigma(2 * i + j)
            Rn = QQ((-1) ** i) * sigma(j)
            d = L - Rn
            if d != 0:
                all_ok = False
            print("(%d,%d) | %5s | %5s | %s" % (i, j, L, Rn, d))

    print("")
    print("(B) (sigma_{2i+j})_{i,j<n} = u w^T ;  (C) det(A+s uw^T) = det(A)(1+s w^T A^{-1} u)")
    print("  n | (B) equal | det(A)=HH^E | (C) dif (should be 0)")
    print("----+-----------+-------------+----------------------")
    for n in range(1, N + 1):
        A = matrix(QQ, n, n, lambda i, j: a[2 * i + j])
        u = vector(QQ, [(-1) ** i for i in range(n)])
        w = vector(QQ, [sigma(j) for j in range(n)])
        Pert = matrix(QQ, n, n, lambda i, j: sigma(2 * i + j))
        outer = matrix(QQ, n, n, lambda i, j: u[i] * w[j])
        B_ok = (Pert == outer)
        detA = A.det()
        detA_ok = (detA == HE(n))
        # (C): det over QQ[s]
        As = matrix(R, n, n, lambda i, j: R(a[2 * i + j]) + s * R(u[i] * w[j]))
        L = As.det()
        Rn = R(detA) * (1 + s * R(w * A.inverse() * u))
        d = L - Rn
        if not (B_ok and detA_ok) or d != 0:
            all_ok = False
        print("%3d | %9s | %11s | %s" % (n, B_ok, detA_ok, d))

    print("")
    if all_ok:
        print("Resume: (A),(B),(C) all hold for indices up to %d. Formula checked." % N)
    else:
        print("Resume: some check failed. Formula NOT checked!")
