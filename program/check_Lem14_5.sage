# check_Lem14_5.sage -- verification of Lemma 14.5 (lem:ds-cauchy) of d08xcosds.tex.
#
# The Cauchy determinant that collapses the structure determinant at s=1.
# For integers M>=1 and N>=0,  Lemma 14.5 (eq:ds-cauchy) asserts
#   det( 1 / (4(M+r-a)^2 - 1) )_{0<=r,a<=N-1}
#     = (-1)^{C(N,2)} 4^{-N}
#       * ( prod_{k=0}^{N-1} k! ) ( prod_{k=1}^{N} k! )
#       / prod_{k=1}^{N} ( M^2 - (2k-1)^2/4 )^{N+1-k}.
#
# M is kept SYMBOLIC (indeterminate over QQ), so each identity is checked as a
# rational-function identity in M, valid for all M.  Indexed by N = 1..Nmax.
#
# Run:  sage check_Lem14_5.sage Nmax

import sys

RM = PolynomialRing(QQ, 'M'); M = RM.gen()
FM = RM.fraction_field()
Mv = FM(M)

def LHS(Nn):
    return matrix(FM, Nn, Nn,
                  lambda r, a: 1 / (4 * (Mv + r - a) ** 2 - 1)).det()

def RHS(Nn):
    return ((-1) ** binomial(Nn, 2) * QQ(4) ** (-Nn)
            * prod(factorial(k) for k in range(0, Nn))
            * prod(factorial(k) for k in range(1, Nn + 1))
            / prod((Mv ** 2 - QQ((2 * k - 1) ** 2) / 4) ** (Nn + 1 - k)
                   for k in range(1, Nn + 1)))

if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("Lemma 14.5 (lem:ds-cauchy) of d08xcosds.tex   [M kept symbolic]")
    print("det( 1/(4(M+r-a)^2-1) )_{0<=r,a<N} = closed form.")
    print("")
    print("  N |  dif = LHS(N) - RHS(N)")
    print("----+-----------------------")

    all_zero = True
    for Nn in range(1, Nmax + 1):
        d = FM(LHS(Nn) - RHS(Nn))
        if d != 0:
            all_zero = False
        print("%3d | %s" % (Nn, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for N = 1..%d. Formula checked." % Nmax)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
