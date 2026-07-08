# check_Thm14_1.sage -- verification of Theorem 14.1 (thm:dblshift) of d08xcosds.tex.
#
# The DOUBLE shift of the secant-number family f(x) = (1+x)/cos(x)^{s+1}, at s=1.
# The coefficient sequence a=(a_n) of f is
#   a_{2k}   = E_{2k}^{(s+1)} = (2k)! [x^{2k}] cos(x)^{-(s+1)},
#   a_{2k+1} = (2k+1) a_{2k}                       (since f = (1+x) cos^{-(s+1)}).
# The double-shifted dilated determinant is
#   HH_n^{(2)} = det(a_{2i+j+2})_{0<=i,j<n},
# with HH_n = det(a_{2i+j}) and HH_n^{(1)} = det(a_{2i+j+1}).
# Theorem 14.1 (eq:ds-target) asserts, at s=1 and for all n>=1,
#   HH_n^{(2)}|_{s=1}
#     = 2^n (n!)^2 HH_n|_{s=1}                                     [form A]
#     = 2^{n+C(n,2)} (n!)^3 ((n-1)!!)^2 prod_{k=1}^{n-2} (k!!)^6    [form B, closed]
#   and equivalently  HH_n^{(2)}/HH_n^{(1)}|_{s=1} = 2^n (n!!)^2    [form C].
#
# The three forms are checked independently (exact integer arithmetic at s=1).
#
# Run:  sage check_Thm14_1.sage N

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()

def dfact(m):
    # double factorial m!! ; dfact(-1)=dfact(0)=1
    r = ZZ(1); k = m
    while k > 1:
        r *= k; k -= 2
    return r

def build_a(M):
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2)
    x = Rx.gen()
    logcos = cos(x).log()
    Ssec = (-(s + 1) * logcos).exp()                 # cos^{-(s+1)}
    E = [factorial(2 * k) * Ssec[2 * k] for k in range(M // 2 + 1)]
    a = []
    for m in range(M + 1):
        if m % 2 == 0:
            a.append(E[m // 2])
        else:
            a.append(m * E[(m - 1) // 2])            # a_{2k+1} = (2k+1) E_{2k}
    return a

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N + 2)                            # max index 2(n-1)+(n-1)+2 = 3n-1
    a1 = [ZZ(Rs(c).subs(s=1)) for c in a]             # specialise s=1

    def HH(shift, n):
        return matrix(ZZ, n, n, lambda i, j: a1[2 * i + j + shift]).det()

    print("Theorem 14.1 (thm:dblshift) of d08xcosds.tex   [double shift, s=1]")
    print("HH_n^{(2)} vs  2^n(n!)^2 HH_n  vs  closed product  vs  2^n(n!!)^2 HH_n^{(1)}")
    print("")
    print("  n |          HH_n^{(2)}|_{s=1}           | difA | difB | difC")
    print("----+--------------------------------------+------+------+-----")

    all_zero = True
    for n in range(1, N + 1):
        L = HH(2, n)
        A = 2 ** n * factorial(n) ** 2 * HH(0, n)
        B = (2 ** (n + binomial(n, 2)) * factorial(n) ** 3 * dfact(n - 1) ** 2
             * prod(dfact(k) ** 6 for k in range(1, n - 1)))
        C = 2 ** n * dfact(n) ** 2 * HH(1, n)
        dA, dB, dC = L - A, L - B, L - C
        if dA != 0 or dB != 0 or dC != 0:
            all_zero = False
        print("%3d | %36s | %4s | %4s | %s" % (n, L, dA, dB, dC))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
