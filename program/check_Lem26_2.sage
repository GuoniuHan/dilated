# check_Lem26_2.sage -- verification of Lemma 26.2 (lem:rhored) of d50chapoton.tex.
#
# Reduction of rho on powers of (1+x).  With A_k = rho(x^k), B_k = rho((1+x)x^k)
# and rho the Chapoton--Han index-lowering evaluation (see check_Thm26_1.sage),
# the lemma asserts, for all i,m >= 0,
#   rho( x^i (1+x)^{2m}   ) = (-2)^m sum_{l=0}^m C(m,l)(-2)^l A_{i+l},
#   rho( x^i (1+x)^{2m+1} ) = (-2)^m sum_{l=0}^m C(m,l)(-2)^l B_{i+l}.
#
# Everything is exact over QQ.  For each pair (i,m) the left side is computed by
# applying rho directly; the right side uses A_{i+l}=rho(x^{i+l}), B similarly.
#
# Indexed by pairs; "sage check_Lem26_2.sage N" checks all 0 <= i,m <= N.
#
# Run:  sage check_Lem26_2.sage N
# prints, for each (i,m), the even/odd difs, then a resume.

import sys

Rx = PolynomialRing(QQ, 'x')
x = Rx.gen()

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

    print("Lemma 26.2 (lem:rhored) of d50chapoton.tex   [exact over QQ]")
    print("rho(x^i(1+x)^{2m})=(-2)^m sum C(m,l)(-2)^l A_{i+l};  odd power -> B")
    print("")
    print("(i,m) | rho(x^i(1+x)^2m) | dif(even) | rho(x^i(1+x)^{2m+1}) | dif(odd)")
    print("------+------------------+-----------+----------------------+---------")

    all_zero = True
    for i in range(N + 1):
        for m in range(N + 1):
            Le = rho(x**i * (1 + x)**(2*m))
            Re = (-2)**m * sum(binomial(m, l) * (-2)**l * A(i + l) for l in range(m + 1))
            Lo = rho(x**i * (1 + x)**(2*m + 1))
            Ro = (-2)**m * sum(binomial(m, l) * (-2)**l * B(i + l) for l in range(m + 1))
            de, do = Le - Re, Lo - Ro
            if de != 0 or do != 0:
                all_zero = False
            print("(%d,%d) | %16s | %9s | %20s | %s" % (i, m, Le, de, Lo, do))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= i,m <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
