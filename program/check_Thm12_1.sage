# check_Thm12_1.sage -- verification of Theorem 12.1 (thm:allstar) of d06xcos.tex.
#
# The secant-number family (1+x)/cos(x)^{s+1} for GENERAL s.  With
#   a_{2k} = E^{(s+1)}_{2k},   a_{2k+1} = (2k+1) E^{(s+1)}_{2k},
# where E^{(s+1)}_{2k} = (2k)! [x^{2k}] sec(x)^{s+1} are the generalised secant
# numbers (polynomials in s), the dilated Hankel determinant
#   HH_n(f) = det(a_{2i+j})_{0<=i,j<=n-1}
# is claimed to equal
#   HH_n(f) = c_n * prod_{i=1}^{n-1} (s+1)_i = c_n * prod_{j=1}^{n-1} (s+j)^{n-j},
#   c_n = 2^C(n,2) ((n-1)!!)^2 prod_{k=1}^{n-2}(k!!)^6 / prod_{i=1}^{n-1} i!.
#
# STRONGEST check: s is kept as an INDETERMINATE, so LHS and RHS are compared as
# polynomials in s over QQ (dif is the zero polynomial).
#
# Run:  sage check_Thm12_1.sage N
# prints, for n = 1..N, whether dif = LHS(n)-RHS(n) is the zero polynomial,
# then a resume stating whether every dif is 0.

import sys

R = PolynomialRing(QQ, 's')
s = R.gen()

def df(k):
    # double factorial k!!  (k>=0); 0!! = (-1)!! = 1
    return prod(range(k, 0, -2))

def E_moments(Kmax):
    # E^{(s+1)}_{2k} = (2k)! [x^{2k}] sec(x)^{s+1}, k = 0..Kmax, as elements of QQ[s].
    # G = sec^{s+1} solves (log G)' = (s+1) tan x, i.e. G' = (s+1) tan(x) G;
    # this gives an exact coefficient recurrence over QQ[s].
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

def poch(base, i):
    # rising factorial (base)_i = base (base+1) ... (base+i-1)
    return prod(base + j for j in range(i))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    E = E_moments(2 * N)                     # index needed: up to 3n-3 < 2*(2N)

    def a(m):
        k = m // 2
        return E[k] if m % 2 == 0 else (2 * k + 1) * E[k]

    def LHS(n):
        return matrix(R, n, n, lambda i, j: a(2 * i + j)).det()

    def cn(n):
        return QQ(2) ** binomial(n, 2) * df(n - 1) ** 2 \
            * prod(df(k) ** 6 for k in range(1, n - 1)) \
            / prod(factorial(i) for i in range(1, n))

    def RHS(n):
        return R(cn(n)) * prod(poch(s + 1, i) for i in range(1, n))

    print("Theorem 12.1 (thm:allstar) of d06xcos.tex   [exact, s indeterminate over QQ]")
    print("E^{(s+1)}_2 = %s,  E^{(s+1)}_4 = %s" % (E[1], E[2]))
    print("")
    print("  n | deg_s | dif = LHS(n) - RHS(n)  (should be 0);  LHS(n) factored")
    print("----+-------+----------------------------------------------------------")

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
