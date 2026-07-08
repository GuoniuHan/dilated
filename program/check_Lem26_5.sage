# check_Lem26_5.sage -- verification of Lemma 26.5 (lem:rhomom) of d50chapoton.tex.
#
# Moment representation of the evaluation rho.  Let S[y^j]=E_{2j} and
# T*[y^j]=(2j+1)E_{2j} be the secant functionals of Section 11, and
# zeta = (x-1)^2/x.  A palindromic polynomial of even index 2e is x^e pi(zeta),
# of odd index 2e+1 is (1+x)x^e sigma(zeta).  With zeta* = -4/(1+y), the lemma asserts
#   rho( x^e pi(zeta) )         = S [ ((1+y)/2)^e pi(zeta*)  ],
#   rho( (1+x) x^e sigma(zeta) )= T*[ ((1+y)/2)^e sigma(zeta*) ],
# and in particular (pi=sigma=1, e=k)
#   A_k = rho(x^k)      = 2^{-k} sum_{j=0}^k C(k,j)   E_{2j},
#   B_k = rho((1+x)x^k) = 2^{-k} sum_{j=0}^k C(k,j)(2j+1)E_{2j}.
#
# Everything is exact over QQ.  Two checks:
#   (A) A_k, B_k via rho vs the E-sum formulas [eq:ABmom];
#   (B) the general representation with pi=sigma=zeta^p (a basis), using the
#       polynomial identity x^e zeta^p = x^{e-p}(x-1)^{2p}, over 0 <= p <= e.
#
# Run:  sage check_Lem26_5.sage N
# prints two tables, then a resume.

import sys

Rx = PolynomialRing(QQ, 'x')
x = Rx.gen()
Ry = PolynomialRing(QQ, 'y')
y = Ry.gen()
Fy = Ry.fraction_field()

def index(P):
    return P.degree() + P.valuation()

def N0(P):
    d = index(P)
    num = (x**d + 1) * P(1) - 2 * P
    q, r = num.quo_rem((x - 1)**2)
    assert r == 0, "N0 division not exact"
    return q

def rho(P):
    if P == 0:
        return QQ(0)
    d = index(P)
    for _ in range(d // 2):
        P = N0(P)
    return P[0]

A = lambda k: rho(x**k)
B = lambda k: rho((1 + x) * x**k)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    # secant numbers and functionals
    Rt = PowerSeriesRing(QQ, 'tt', default_prec=2 * N + 3)
    tt = Rt.gen()
    sec = 1 / cos(tt)
    E = lambda j: factorial(2 * j) * sec[2 * j]
    Smom = lambda j: E(j)
    Tsmom = lambda j: (2 * j + 1) * E(j)
    applyf = lambda poly, mom: sum(c * mom(j) for j, c in enumerate(poly.list()))
    zeta_star = -4 / (1 + y)

    print("Lemma 26.5 (lem:rhomom) of d50chapoton.tex   [exact over QQ]")
    print("")

    all_zero = True

    print("(A) A_k, B_k  vs  2^{-k} sum C(k,j) E_{2j} (resp. (2j+1)E_{2j})")
    print("  k |   A_k   | dif |   B_k   | dif")
    print("----+---------+-----+---------+-----")
    for k in range(N + 1):
        ra = 2**(-k) * sum(binomial(k, j) * E(j) for j in range(k + 1))
        rb = 2**(-k) * sum(binomial(k, j) * (2*j + 1) * E(j) for j in range(k + 1))
        da, db = A(k) - ra, B(k) - rb
        if da != 0 or db != 0:
            all_zero = False
        print("%3d | %7s | %3s | %7s | %s" % (k, A(k), da, B(k), db))

    print("")
    print("(B) general rho(x^e pi(zeta)) = S[((1+y)/2)^e pi(zeta*)], pi=zeta^p (and T* odd)")
    print("(e,p) | even dif | odd dif")
    print("------+----------+--------")
    for e in range(N + 1):
        for p in range(e + 1):
            # x^e zeta^p = x^{e-p}(x-1)^{2p}
            Pe = x**(e - p) * (x - 1)**(2*p)
            Po = (1 + x) * x**(e - p) * (x - 1)**(2*p)
            rhs = Fy(((1 + y)/2)**e * zeta_star**p)
            assert rhs.denominator() == 1, "RHS not a polynomial in y"
            Re = applyf(rhs.numerator(), Smom)
            Ro = applyf(rhs.numerator(), Tsmom)
            de, do = rho(Pe) - Re, rho(Po) - Ro
            if de != 0 or do != 0:
                all_zero = False
            print("(%d,%d) | %8s | %s" % (e, p, de, do))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for indices <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
