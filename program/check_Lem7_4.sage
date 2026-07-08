# check_Lem7_4.sage -- verification of Lemma 7.4 (lem:bindetSigma) of d04euler.tex.
#
# Binomial determinant, general parameter.  For every N >= 1, q in {N, N+1}, and
# arbitrary a,
#   det( binomial(a, q + r - m) )_{0<=r,m<=N-1}
#     = prod_{r=1}^{q} prod_{c=1}^{N} (a + c - r) / ((q - r) + (N - c) + 1).
#
# The parameter a is kept SYMBOLIC (an indeterminate over QQ), so each identity
# is a polynomial identity in a, valid for all a.  For matrix size N we check
# BOTH admissible column offsets q = N and q = N+1.
#
# Indexed by pairs (N, q); "sage check_Lem7_4.sage NMAX" checks N = 1..NMAX.
#
# Run:  sage check_Lem7_4.sage NMAX
# prints, for each (N, q), LHS, RHS and dif = LHS-RHS,
# then a resume stating whether every dif is 0.

import sys

Ra = PolynomialRing(QQ, 'a')
a = Ra.gen()

def binom_a(k):
    # binomial(a, k) = a(a-1)...(a-k+1)/k!  for integer k >= 0;  0 for k < 0
    if k < 0:
        return Ra(0)
    return prod(a - j for j in range(k)) / factorial(k)

def LHS(Nsize, q):
    return matrix(Ra, Nsize, Nsize,
                  lambda r, m: binom_a(q + r - m)).det()

def RHS(Nsize, q):
    return prod((a + c - r) / ((q - r) + (Nsize - c) + 1)
                for r in range(1, q + 1) for c in range(1, Nsize + 1))

if __name__ == "__main__":
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("Lemma 7.4 (lem:bindetSigma) of d04euler.tex   [a kept symbolic]")
    print("")
    print("(N,q) |                     LHS                     |  dif")
    print("------+---------------------------------------------+-----")

    all_zero = True
    for Nsize in range(1, NMAX + 1):
        for q in (Nsize, Nsize + 1):
            L = LHS(Nsize, q)
            R = RHS(Nsize, q)
            d = Ra(L - R)
            if d != 0:
                all_zero = False
            print("(%d,%d) | %43s | %s" % (Nsize, q, L, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for N = 1..%d (q in {N,N+1}). Formula checked." % NMAX)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
