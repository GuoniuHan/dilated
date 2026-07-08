# check_Lem10_10.sage -- verification of Lemma 10.10 (lem:t3trace) of d05eulershift.tex.
#
# Let Mbar = ( C(nbar+r, m) )_{0<=r,m<=nund-1}, Dm = diag(0,1,...,nund-1), and
# Da = diag(a_0,...,a_{nund-1}) an ARBITRARY diagonal matrix.  Then (eq:t3traceid)
#   tr( Mbar^{-1} Da Mbar Dm ) = sum_{r=0}^{nund-1} a_r (nbar+r) - nbar a_0 nund.
#
# The a_r are kept as indeterminates, so each n gives a polynomial identity in
# (a_0,...,a_{nund-1}).
#
# Run:  sage check_Lem10_10.sage N

import sys

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    print("Lemma 10.10 (lem:t3trace) of d05eulershift.tex   [a_r symbolic]")
    print("")
    print("  n | nund |            dif = LHS - RHS (identity in a_0..a_{nund-1})")
    print("----+------+---------------------------------------------------------")

    all_zero = True
    for n in range(1, N + 1):
        nbar = (n + 1)//2; nund = n//2
        if nund == 0:
            print("%3d | %4d | 0   (empty)" % (n, nund))
            continue
        A = PolynomialRing(QQ, ['a%d' % i for i in range(nund)])
        avars = A.gens()
        Mbar = matrix(QQ, nund, nund, lambda r, m: binomial(nbar + r, m))
        Minv = Mbar.inverse()
        Da = diagonal_matrix(A, list(avars))
        Dm = diagonal_matrix(A, [A(i) for i in range(nund)])
        MbarA = Mbar.change_ring(A); MinvA = Minv.change_ring(A)
        LHS = (MinvA * Da * MbarA * Dm).trace()
        RHS = sum(avars[r]*(nbar + r) for r in range(nund)) - nbar*avars[0]*nund
        d = A(LHS - RHS)
        if d != 0: all_zero = False
        print("%3d | %4d | %s" % (n, nund, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Identity checked." % N)
    else:
        print("Resume: some dif != 0. Identity NOT checked!")
