# check_Lem12_3.sage -- verification of Lemma 12.3 (lem:rankdrop) of d06xcos.tex.
#
# Rank drop at negative integers for f(x) = (1+x)/cos(x)^{s+1}.  The lemma asserts:
# for each integer p with 1 <= p <= n-1, the matrix (a_{2i+j})_{0<=i,j<n} has
# rank AT MOST p at s = -p; consequently
#   (s+p)^{n-p} | HH_n(f),    and    prod_{j=1}^{n-1} (s+j)^{n-j} | HH_n(f).
#
# Two exact checks:
#   (A) numeric rank of the matrix at s = -p (over QQ) is <= p;
#   (B) HH_n(f) (as a polynomial in s over QQ) is divisible by (s+p)^{n-p},
#       and by the full product prod_{j=1}^{n-1}(s+j)^{n-j} (of degree C(n,2)).
#
# At s = -p the even weight is cos(x)^{p-1}, so a_{2k} = (2k)![x^{2k}] cos(x)^{p-1}
# and a_{2k+1} = (2k+1) a_{2k}; these are exact rationals over QQ.
#
# Run:  sage check_Lem12_3.sage N
# prints, for each (n,p), the rank and the divisibility flags, then a resume.

import sys

R = PolynomialRing(QQ, 's')
s = R.gen()

def E_moments(Kmax):
    # generalised secant numbers E^{(s+1)}_{2k} over QQ[s]  (for the polynomial HH_n)
    Nser = 2 * Kmax + 2
    St = PowerSeriesRing(QQ, 'x', default_prec=Nser + 3)
    x = St.gen()
    tanser = sin(x) / cos(x)
    t = [QQ(tanser[i]) for i in range(Nser + 1)]
    g = [R(0)] * (Nser + 1)
    g[0] = R(1)
    for j in range(Nser):
        acc = R(0)
        for i in range(1, j + 1):
            acc += t[i] * g[j - i]
        g[j + 1] = (s + 1) * acc / (j + 1)
    return [factorial(2 * k) * g[2 * k] for k in range(Kmax + 1)]

def a_at_negp(Kmax, p):
    # a_{2k}, a_{2k+1} at s = -p over QQ: even weight cos(x)^{p-1}
    Nser = 2 * Kmax + 2
    St = PowerSeriesRing(QQ, 'x', default_prec=Nser + 3)
    x = St.gen()
    w = cos(x) ** (p - 1)
    Eneg = [factorial(2 * k) * QQ(w[2 * k]) for k in range(Kmax + 1)]
    def a(m):
        k = m // 2
        return Eneg[k] if m % 2 == 0 else (2 * k + 1) * Eneg[k]
    return a

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    E = E_moments(2 * N)

    def a_sym(m):
        k = m // 2
        return E[k] if m % 2 == 0 else (2 * k + 1) * E[k]

    def HH_poly(n):
        return matrix(R, n, n, lambda i, j: a_sym(2 * i + j)).det()

    print("Lemma 12.3 (lem:rankdrop) of d06xcos.tex   [exact over QQ]")
    print("")
    print("(A) rank of (a_{2i+j})_{0<=i,j<n} at s=-p  should be <= p")
    print("(B) (s+p)^{n-p} | HH_n   and   prod_j (s+j)^{n-j} | HH_n")
    print("")
    print("(n,p) | rank | p | rank<=p | (s+p)^{n-p}|HH | full-prod|HH")
    print("------+------+---+---------+---------------+-------------")

    all_ok = True
    for n in range(1, N + 1):
        H = HH_poly(n)
        fullprod = prod((s + j) ** (n - j) for j in range(1, n))
        full_div = (H % fullprod == 0)
        if not full_div:
            all_ok = False
        for p in range(1, n):
            a = a_at_negp(2 * n, p)
            A = matrix(QQ, n, n, lambda i, j: a(2 * i + j))
            rk = A.rank()
            rank_ok = (rk <= p)
            div = ((s + p) ** (n - p)).divides(H)
            if not (rank_ok and div):
                all_ok = False
            print("(%d,%d) | %4d | %d | %7s | %13s | %s"
                  % (n, p, rk, p, rank_ok, div, full_div))

    print("")
    if all_ok:
        print("Resume: rank<=p and all divisibilities hold for n=1..%d. Formula checked." % N)
    else:
        print("Resume: some check failed. Formula NOT checked!")
