# check_Lem5_3.sage -- verification of Lemma 5.3 (lem:mixed) of d03invo.tex.
#
# Mixed Hermite moments.  Let L be the Gaussian moment functional
#   L[z^m] = (m-1)!! if m even, 0 if m odd,
# and He_m the monic (probabilist) Hermite polynomials.  Then for all p,q >= 0
#   L[ He_p(c+z) He_q(z) ] = q! binomial(p,q) c^{p-q}
# (in particular 0 for q > p).
#
# c is kept SYMBOLIC, so each identity is checked as a polynomial identity in c.
# Indexed by the pair (p,q); "sage check_Lem5_3.sage N" checks all 0<=p,q<=N.
#
# Run:  sage check_Lem5_3.sage N
# prints, for each (p,q), LHS, RHS and dif = LHS-RHS,
# then a resume stating whether every dif is 0.

import sys

Rc = PolynomialRing(QQ, 'c'); c = Rc.gen()
Rz = PolynomialRing(Rc, 'z'); z = Rz.gen()

def He(m):
    # monic Hermite: He_0=1, He_1=z, He_{m+1}=z He_m - m He_{m-1}
    h0, h1 = Rz(1), z
    if m == 0: return h0
    for k in range(1, m):
        h0, h1 = h1, z*h1 - k*h0
    return h1

def Lmom(m):
    # L[z^m]
    return Rc(prod(range(m - 1, 0, -2))) if m % 2 == 0 else Rc(0)

def Lfun(P):
    # apply L to a polynomial P in z (coefficients in Rc)
    return sum(coef * Lmom(m) for m, coef in enumerate(P.list()))

def LHS(p, q):
    return Lfun(He(p)(c + z) * He(q)(z))

def RHS(p, q):
    return factorial(q) * binomial(p, q) * c^(p - q) if q <= p else Rc(0)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("Lemma 5.3 (lem:mixed) of d03invo.tex   [c kept symbolic]")
    print("checking all pairs 0 <= p,q <= %d" % N)
    print("")
    print("(p,q) |            LHS            |            RHS            |  dif")
    print("------+--------------------------+--------------------------+-----")

    all_zero = True
    for p in range(N + 1):
        for q in range(N + 1):
            L = LHS(p, q)
            R = RHS(p, q)
            d = L - R
            if d != 0:
                all_zero = False
            print("(%d,%d) | %24s | %24s | %s" % (p, q, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= p,q <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
