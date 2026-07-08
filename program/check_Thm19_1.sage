# check_Thm19_1.sage -- verification of Theorem 19.1 (thm:ell) of d15elliptic.tex.
#
# An elliptic deformation of the Euler numbers.  With the Jacobi elliptic
# functions sn(x,m), cn(x,m) (modulus parameter m), the Euler generating
# function (1+sin x)/cos x is deformed into
#   g_m(x) = (1 + sn(x,m)) / cn(x,m) = sum_{n>=0} E_n(m) x^n/n!,
# whose coefficients a_k = E_k(m) are polynomials in m (E_2=1, E_3=2-m,
# E_4=5-4m, E_5=16-16m+m^2, ...).  Theorem 19.1 asserts, for all m and n>=1,
#   HH_n(g_m) = det(a_{2i+j})_{0<=i,j<n} = ((1-m)/2)^C(n,2) * prod_{k=1}^{n-1} k! (2k)!.
#
# The sn, cn, dn Taylor series are generated over QQ[m] from their defining ODEs
#   sn' = cn dn,  cn' = -sn dn,  dn' = -m sn cn,   sn(0)=0, cn(0)=1, dn(0)=1,
# so m is kept SYMBOLIC (indeterminate over QQ) and each LHS(n)=RHS(n) is a
# polynomial identity in m, valid for ALL m.
#
# Run:  sage check_Thm19_1.sage N
# prints, for n = 1..N, LHS(n), the difference dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

Rm = PolynomialRing(QQ, 'm')
m = Rm.gen()

def jacobi_series(M):
    # Taylor coefficients of sn, cn, dn up to x^M over QQ[m], from the ODEs.
    s = [Rm(0)]; c = [Rm(1)]; d = [Rm(1)]
    for k in range(M):
        conv_cd = sum(c[j] * d[k - j] for j in range(k + 1))   # [x^k](cn*dn)
        conv_sd = sum(s[j] * d[k - j] for j in range(k + 1))   # [x^k](sn*dn)
        conv_sc = sum(s[j] * c[k - j] for j in range(k + 1))   # [x^k](sn*cn)
        s.append(conv_cd / (k + 1))
        c.append(-conv_sd / (k + 1))
        d.append(-m * conv_sc / (k + 1))
    return s, c, d

def build_a(M):
    # a_k = k! [x^k] (1+sn)/cn ,  k = 0..M
    s, c, d = jacobi_series(M)
    Rx = PowerSeriesRing(Rm, 'x', default_prec=M + 1)
    x = Rx.gen()
    sn = sum(s[k] * x ** k for k in range(M + 1))
    cn = sum(c[k] * x ** k for k in range(M + 1))
    g = (1 + sn) / cn
    return [factorial(k) * g[k] for k in range(M + 1)]

def LHS(n, a):
    return matrix(Rm, n, n, lambda i, j: a[2 * i + j]).det()

def RHS(n):
    return ((1 - m) / 2) ** binomial(n, 2) \
        * prod(factorial(k) * factorial(2 * k) for k in range(1, n))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)                         # max index 2(n-1)+(n-1) = 3n-3

    print("Theorem 19.1 (thm:ell) of d15elliptic.tex   [m kept symbolic]")
    print("HH_n(g_m) = det(a_{2i+j}) = ((1-m)/2)^C(n,2) prod_{k=1}^{n-1} k!(2k)!")
    print("")
    print("  n |                       LHS(n) = HH_n(g_m)                       |  dif")
    print("----+----------------------------------------------------------------+-----")

    all_zero = True
    for n in range(1, N + 1):
        L = LHS(n, a)
        R = RHS(n)
        d = Rm(L - R)
        if d != 0:
            all_zero = False
        print("%3d | %62s | %s" % (n, L, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
