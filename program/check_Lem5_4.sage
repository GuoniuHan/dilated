# check_Lem5_4.sage -- verification of Lemma 5.4 (lem:tri) of d03invo.tex.
#
# Triangularising family.  With L the Gaussian moment functional, He_m the monic
# Hermite polynomials, and
#   Phi_i(z) = sum_{k=0}^i binomial(i,k) (-c^2)^k He_{2i-2k}(c+z),
# one has for all i,j >= 0
#   L[ Phi_i(z) He_j(z) ] = j! binomial(i, j-i) (2c)^{2i-j}
# (in particular 0 for j < i, and i! (2c)^i for j = i).
#
# c is kept SYMBOLIC, so each identity is checked as a polynomial identity in c.
# Indexed by the pair (i,j); "sage check_Lem5_4.sage N" checks all 0<=i,j<=N.
#
# Run:  sage check_Lem5_4.sage N
# prints, for each (i,j), LHS, RHS and dif = LHS-RHS,
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
    return Rc(prod(range(m - 1, 0, -2))) if m % 2 == 0 else Rc(0)

def Lfun(P):
    return sum(coef * Lmom(m) for m, coef in enumerate(P.list()))

def Phi(i):
    return sum(binomial(i, k) * (-c^2)^k * He(2*i - 2*k)(c + z)
               for k in range(i + 1))

def LHS(i, j):
    return Lfun(Phi(i) * He(j)(z))

def RHS(i, j):
    return factorial(j) * binomial(i, j - i) * (2*c)^(2*i - j) if i <= j <= 2*i else Rc(0)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("Lemma 5.4 (lem:tri) of d03invo.tex   [c kept symbolic]")
    print("checking all pairs 0 <= i,j <= %d" % N)
    print("")
    print("(i,j) |            LHS            |            RHS            |  dif")
    print("------+--------------------------+--------------------------+-----")

    all_zero = True
    for i in range(N + 1):
        for j in range(N + 1):
            L = LHS(i, j)
            R = RHS(i, j)
            d = L - R
            if d != 0:
                all_zero = False
            print("(%d,%d) | %24s | %24s | %s" % (i, j, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= i,j <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
