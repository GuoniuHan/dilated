# check_Prop15_1.sage -- verification of Proposition 15.1 (prop:sin3) of d10runkone.tex.
#
# Rank-one perturbation of the shifted Euler family:
#   f_s(x) = (sin x + 1)/cos^2 x + s sin x.
# Its coefficient sequence is tilde a_n = a_n + s sigma_n, where
#   a_n = n! [x^n] (sin x + 1)/cos^2 x = E_{n+1}   (shifted Euler / zigzag),
#   sigma_{2k} = 0,  sigma_{2k+1} = (-1)^k        (coefficients of sin x).
# The proposition claims, for all s and all n >= 1,
#   HH_n(f_s) = (1 - s C(n,2)) * HH_n^E,   HH_n^E = prod_{k=1}^{n-1} (k!)^2 (2k+1)!!.
#
# STRONGEST check: s is kept as an INDETERMINATE; LHS = det(tilde a_{2i+j}) and RHS
# are compared as polynomials in s over QQ (dif is the zero polynomial).
#
# Run:  sage check_Prop15_1.sage N
# prints, for n=1..N, dif = LHS(n)-RHS(n) (should be 0) and LHS factored, then a resume.

import sys

R = PolynomialRing(QQ, 's')
s = R.gen()

def df(k):
    return prod(range(k, 0, -2))          # double factorial k!!

def build_a(M):
    # a_n = n! [x^n] (sin x + 1)/cos^2 x, exact over QQ, n = 0..M
    Rx = PowerSeriesRing(QQ, 'x', default_prec=M + 2)
    x = Rx.gen()
    g = (sin(x) + 1) / cos(x) ** 2
    return [factorial(m) * g[m] for m in range(M + 1)]

def sigma(m):
    # coefficients of sin x: sigma_{2k}=0, sigma_{2k+1}=(-1)^k
    if m % 2 == 0:
        return QQ(0)
    return QQ((-1) ** ((m - 1) // 2))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)                     # max index 3n-3

    def ta(m):
        return R(a[m]) + s * R(sigma(m))

    def LHS(n):
        return matrix(R, n, n, lambda i, j: ta(2 * i + j)).det()

    def HE(n):
        return prod(factorial(k) ** 2 * df(2 * k + 1) for k in range(1, n))

    def RHS(n):
        return (1 - s * binomial(n, 2)) * R(HE(n))

    print("Proposition 15.1 (prop:sin3) of d10runkone.tex   [exact, s indeterminate over QQ]")
    print("a = %s ...  (a_n = E_{n+1})" % a[:8])
    print("")
    print("  n | deg_s | dif = LHS(n)-RHS(n) (should be 0);  LHS(n) factored")
    print("----+-------+-------------------------------------------------------")

    all_zero = True
    for n in range(1, N + 1):
        L = LHS(n)
        Rn = RHS(n)
        d = L - Rn
        if d != 0:
            all_zero = False
        fac = L.factor() if L != 0 else "0"
        print("%3d | %5s | dif=%s ;  LHS=%s" % (n, L.degree(), d, fac))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
