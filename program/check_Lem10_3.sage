# check_Lem10_3.sage -- verification of Lemma 10.3 (lem:t3cauchybinet) of d05eulershift.tex.
#
# At t=3, with kt_{i,m} the mixed connection of Lemma 10.1, b=(s-3)/2, and
#   N_{r,m} = (2m+1) C(b,d) + (s+2m+2) C(b,d-1),   d = nbar+r-m,   0<=r,m<=nund-1,
# the lemma asserts three things, checked here for each n:
#   (A) eq:t3rowcol:  det(kt_{nbar+r,m}) = (prod_r (2nbar+2r)!)(prod_m 1/(2m+1)!) det N
#   (B) factorisation N = V B, with V the nund x (nund+1) Pascal block
#         V_{r,k}=C(b, nbar+r-k) and B the (nund+1) x nund bidiagonal matrix
#         B_{k,m}=(2m+1)[k=m]+(s+2m+2)[k=m+1]
#   (C) eq:t3cb (Cauchy-Binet):  det N = sum_{p=0}^{nund} c_p D_p
#         c_p = (prod_{m<p}(2m+1))(prod_{m>=p}(s+2m+2)),  D_p = det V with col p deleted.
#
# s symbolic (t=3); every dif = 0 is a polynomial identity in s.
#
# Run:  sage check_Lem10_3.sage N

import sys

R = PolynomialRing(QQ, 's'); s = R.gen(); F = R.fraction_field()

def gbin(a, d): return F(0) if d < 0 else prod(a - i for i in range(d)) / factorial(d)

def kt_formula(i, m):
    d = i - m
    return (factorial(2*i)/factorial(2*m) * gbin((s-1)/2, d)
            + (s+1)*factorial(2*i)/factorial(2*m+1) * gbin((s-3)/2, d-1))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    b = (s - 3)/2

    print("Lemma 10.3 (lem:t3cauchybinet) of d05eulershift.tex   [s symbolic, t=3]")
    print("Three checks per n:  (A) row/col extraction,  (B) N=VB,  (C) Cauchy-Binet.")
    print("")
    print("  n | difA=det(kt)-fac*detN | difB=max|N-VB| | difC=detN-sum c_p D_p")
    print("----+-----------------------+----------------+----------------------")

    all_zero = True
    for n in range(1, N + 1):
        nbar = (n + 1)//2; nund = n//2
        if nund == 0:
            print("%3d | %21s | %14s | %s   (empty: nund=0)" % (n, 0, 0, 0))
            continue

        # kt determinant and row/col factor
        ktdet = matrix(F, nund, nund, lambda r, m: kt_formula(nbar+r, m)).det()
        rowcol = (prod(factorial(2*nbar+2*r) for r in range(nund))
                  * prod(1/factorial(2*m+1) for m in range(nund)))
        Nmat = matrix(F, nund, nund,
                      lambda r, m: (2*m+1)*gbin(b, nbar+r-m) + (s+2*m+2)*gbin(b, nbar+r-m-1))
        difA = F(ktdet - rowcol*Nmat.det())

        # N = V B
        V = matrix(F, nund, nund+1, lambda r, k: gbin(b, nbar+r-k))
        B = matrix(F, nund+1, nund,
                   lambda k, m: (2*m+1) if k == m else ((s+2*m+2) if k == m+1 else 0))
        difB_zero = all(F(e) == 0 for e in (Nmat - V*B).list())

        # Cauchy-Binet sum
        tot = F(0)
        for p in range(nund+1):
            cp = prod(2*m+1 for m in range(p)) * prod(s+2*m+2 for m in range(p, nund))
            cols = [k for k in range(nund+1) if k != p]
            Dp = matrix(F, nund, nund, lambda r, jj: gbin(b, nbar+r-cols[jj])).det()
            tot += cp*Dp
        difC = F(Nmat.det() - tot)

        if difA != 0 or not difB_zero or difC != 0: all_zero = False
        print("%3d | %21s | %14s | %s" % (n, difA, "0" if difB_zero else "NONZERO", difC))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
