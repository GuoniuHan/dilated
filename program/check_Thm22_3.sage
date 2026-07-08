# check_Thm22_3.sage -- verification of Theorem 22.3 (thm:algsqshift) of d19algesq.tex.
#
# All shifts of the squared algebraic family f(x) = (1+x)^2 (1-x^2)^{-s/2}, with
#   b_{2k}   = mu_k + 2k(2k-1) mu_{k-1},   b_{2k+1} = 2(2k+1) mu_k,
#   mu_k = (2k)!/k! (s/2)_k,   f_1 = (1+x)(1-x^2)^{-s/2}  (moments a_m, Prop 22.1).
# Let HH_n^{(r)}(f) = det(b_{2i+j+r}), HH_n^{(r)}(f_1) = det(a_{2i+j+r}), rho=floor(r/2).
# The theorem asserts:
#   (i)  r odd:   HH_n^{(r)}(f) = 2^n HH_n^{(r)}(f_1)
#                 = 2^{C(n,2)+n} prod_{k=1}^{n-1} k! prod_{i=0}^{n-1} (2i+r) mu_{i+rho}(s).
#   (ii) r=2rho even (s!=3):  HH_n^{(r)}(f) = 2^{n-1} HH_n^{(r)}(f_1) * P_{n,rho}/Theta_{n,rho},
#        Theta_{n,rho} = prod_{i=0}^{n-1}(s-2+2rho+2i),
#        P_{n,rho} = [ (s-2)(2n+2rho-1)!!/(2rho-1)!! + (s-4) Theta_{n,rho} ] / (s-3),
#        a monic degree-n polynomial; and for rho>=1 the fully explicit form
#        HH_n^{(2rho)}(f) = 2^{C(n,2)-1} prod k! prod (2i+2rho)!/(i+rho)! (s/2)_{i+rho-1} * P_{n,rho}.
#
# All checks are exact polynomial identities in the indeterminate s over QQ.
# For (ii) the P/Theta form is verified in cross-multiplied (polynomial) form
# LHS*Theta = 2^{n-1} HH^{(r)}(f_1) P, valid for every s.
#
# Run:  sage check_Thm22_3.sage N   (odd r in {1,3}, even r in {2,4})

import sys

R = PolynomialRing(QQ, 's')
s = R.gen()

def df(k):
    return prod(range(k, 0, -2))          # double factorial k!! (df(-1)=1)

def poch(base, i):
    return prod(base + j for j in range(i))

def mu(k):
    if k < 0:
        return R(0)
    return R(factorial(2 * k) / factorial(k)) * poch(s / 2, k)

def a(m):
    k = m // 2
    return mu(k) if m % 2 == 0 else m * mu(k)

def b(m):
    k = m // 2
    if m % 2 == 0:
        return mu(k) + (2 * k) * (2 * k - 1) * mu(k - 1)
    return 2 * m * mu(k)                   # 2(2k+1)mu_k, m=2k+1

def Theta(n, rho):
    return prod(s - 2 + 2 * rho + 2 * i for i in range(n))

def Ppoly(n, rho):
    DF = df(2 * n + 2 * rho - 1) / df(2 * rho - 1)
    numer = (s - 2) * R(DF) + (s - 4) * Theta(n, rho)
    P, rem = numer.quo_rem(s - 3)
    assert rem == 0, "numerator not divisible by (s-3)!"
    return P

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    def HHb(n, r):
        return matrix(R, n, n, lambda i, j: b(2 * i + j + r)).det()

    def HHa(n, r):
        return matrix(R, n, n, lambda i, j: a(2 * i + j + r)).det()

    all_zero = True

    print("Theorem 22.3 (thm:algsqshift) of d19algesq.tex   [exact, s indeterminate]")
    print("")
    print("(i) r ODD:  HH^{(r)}(f) = 2^n HH^{(r)}(f_1) = 2^{C(n,2)+n} prod k! prod (2i+r)mu_{i+rho}")
    print("(n,r) | dif(vs 2^n HH_f1) | dif(vs closed form)")
    print("------+-------------------+--------------------")
    for n in range(1, N + 1):
        for r in [1, 3]:
            rho = r // 2
            L = HHb(n, r)
            d1 = L - QQ(2) ** n * HHa(n, r)
            closed = (QQ(2) ** (binomial(n, 2) + n)
                      * prod(factorial(k) for k in range(1, n))
                      * prod((2 * i + r) * mu(i + rho) for i in range(n)))
            d2 = L - closed
            if d1 != 0 or d2 != 0:
                all_zero = False
            print("(%d,%d) | %17s | %s" % (n, r, d1, d2))

    print("")
    print("(ii) r=2rho EVEN:  LHS*Theta = 2^{n-1} HH^{(r)}(f_1) P   AND explicit (rho>=1) form")
    print("(n,r) | dif(P/Theta, cross-mult) | dif(explicit) | P monic deg n")
    print("------+--------------------------+---------------+--------------")
    for n in range(1, N + 1):
        for r in [2, 4]:
            rho = r // 2
            L = HHb(n, r)
            P = Ppoly(n, rho)
            Th = Theta(n, rho)
            # compact form, cross-multiplied to stay polynomial:
            dC = L * Th - QQ(2) ** (n - 1) * HHa(n, r) * P
            # explicit rho>=1 form:
            expl = (QQ(2) ** (binomial(n, 2) - 1)
                    * prod(factorial(k) for k in range(1, n))
                    * prod(R(factorial(2 * i + 2 * rho) / factorial(i + rho))
                           * poch(s / 2, i + rho - 1) for i in range(n))
                    * P)
            dE = L - expl
            Pok = (P.degree() == n and P.leading_coefficient() == 1)
            if dC != 0 or dE != 0 or not Pok:
                all_zero = False
            print("(%d,%d) | %24s | %13s | %s" % (n, r, dC, dE, Pok))

    print("")
    if all_zero:
        print("Resume: all dif = 0 (i: r=1,3; ii: r=2,4) for n=1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
