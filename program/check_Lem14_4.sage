# check_Lem14_4.sage -- verification of Lemma 14.4 (lem:ds-Delta) of d08xcosds.tex.
#
# Structure coefficients of the classical Christoffel functional S1.  With
# Phat_m the monic S1-orthogonal polynomials (data in Lemma 14.2) and
#   T_i := y Phat_i' - i Phat_i = sum_{a=0}^{i-1} Delta_{i,a} Phat_a,
# Lemma 14.4 (eq:ds-Delta) gives the closed form, for 0 <= a <= i-1,
#   Delta_{i,a} = (-1)^{i-a-1} * 2^{i-a-1}/(2(i-a)-1) * i!/a!
#                 * (2i+1)!!/(2a+1)!! * ( s + (2(i+a)+3)/(2(i-a)+1) ),
# in particular deg_s Delta_{i,a} = 1.
#
# LHS(i,a) is computed directly from the definition (expanding T_i in the
# Phat basis); RHS(i,a) is the closed form.  s is kept SYMBOLIC, so each
# identity is a polynomial identity in s.  Indexed by the pairs (i,a) with
# 1<=i<=N, 0<=a<i.
#
# Run:  sage check_Lem14_4.sage N

import sys

Rs = PolynomialRing(QQ, 's'); s = Rs.gen()
Ry = PolynomialRing(Rs, 'y'); y = Ry.gen()

def dfact(m):
    r = QQ(1); k = m
    while k > 1:
        r *= k; k -= 2
    return r

def cS1(m):   return (4 * m + 3) * s + 8 * m ** 2 + 12 * m + 5
def lamS1(m): return 2 * m * (2 * m + 1) * (s + 2 * m) * (s + 2 * m + 1)

def Phat(mmax):
    P = [Ry(1)]
    if mmax >= 1: P.append(y - cS1(0))
    for m in range(1, mmax):
        P.append((y - cS1(m)) * P[m] - lamS1(m) * P[m - 1])
    return P

def Delta_def(i, P):
    q = y * P[i].derivative() - i * P[i]
    coeffs = {}
    while q != 0 and q.degree() >= 0:
        d = q.degree(); lc = q.leading_coefficient()
        coeffs[d] = lc; q = q - lc * P[d]
    return coeffs

def Delta_closed(i, a):
    return ((-1) ** (i - a - 1) * QQ(2) ** (i - a - 1) / (2 * (i - a) - 1)
            * factorial(i) / factorial(a) * dfact(2 * i + 1) / dfact(2 * a + 1)
            * (s + QQ(2 * (i + a) + 3) / (2 * (i - a) + 1)))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    P = Phat(N + 1)

    print("Lemma 14.4 (lem:ds-Delta) of d08xcosds.tex   [s kept symbolic]")
    print("Delta_{i,a}: definition (expand T_i in Phat basis) vs closed form.")
    print("")
    print(" (i,a) |            LHS = Delta_{i,a}            |  dif")
    print("-------+----------------------------------------+-----")

    all_zero = True
    for i in range(1, N + 1):
        cdict = Delta_def(i, P)
        for a in range(i):
            L = cdict.get(a, Rs(0))
            R = Delta_closed(i, a)
            d = Rs(L - R)
            if d != 0: all_zero = False
            print("(%d,%d) | %38s | %s" % (i, a, L, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for pairs (i,a), i<=%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
