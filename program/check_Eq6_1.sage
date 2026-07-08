# check_Eq6_1.sage -- verification of identity (6.1) (eq:dergauss) of d03invo.tex.
#
# With L the Gaussian moment functional, u = c + z, He_m the monic Hermite
# polynomials, Phi_s the triangularising family of Lemma 5.4, and
#   Psi_i = sum_{s=0}^i (i!/s!) (-2)^{i-s} Phi_s,
# the single-shift triangularity identity is
#   L[ Psi_i * u * He_j(z) ] = 0            for j < i,
#                            = c i! (2c)^i  for j = i.
# (These are the only cases the identity claims; they give the diagonal that
# yields HH_n(f') = prod c i! (2c)^i.)  So we check the triangular set j <= i.
#
# c is kept SYMBOLIC, so each identity is a polynomial identity in c.
# "sage check_Eq6_1.sage N" checks all pairs 0 <= j <= i <= N.
#
# Run:  sage check_Eq6_1.sage N
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

def Phi(s):
    return sum(binomial(s, k) * (-c^2)^k * He(2*s - 2*k)(c + z)
               for k in range(s + 1))

def Psi(i):
    return sum(Integer(factorial(i)) / factorial(s) * (-2)^(i - s) * Phi(s)
               for s in range(i + 1))

def LHS(i, j):
    return Lfun(Psi(i) * (c + z) * He(j)(z))

def RHS(i, j):
    if j < i:  return Rc(0)
    if j == i: return c * factorial(i) * (2*c)^i
    return None  # j > i: not claimed

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("Identity (6.1) (eq:dergauss) of d03invo.tex   [c kept symbolic]")
    print("checking all pairs 0 <= j <= i <= %d (the cases the identity claims)" % N)
    print("")
    print("(i,j) |            LHS            |            RHS            |  dif")
    print("------+--------------------------+--------------------------+-----")

    all_zero = True
    for i in range(N + 1):
        for j in range(i + 1):
            L = LHS(i, j)
            R = RHS(i, j)
            d = L - R
            if d != 0:
                all_zero = False
            print("(%d,%d) | %24s | %24s | %s" % (i, j, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= j <= i <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
