# check_Lem10_8.sage -- verification of Lemma 10.8 (lem:t3rank) of d05eulershift.tex.
#
# At s=-(2k+1) the matrix N of (eq:t3cb),
#   N_{r,m} = (2m+1) C(b, d) + (s+2m+2) C(b, d-1),  b=(s-3)/2,  d=nbar+r-m,
# has rank <= k+1 (for 0<=k<=nund-1); equivalently corank >= nund-1-k, so
#   (s+2k+1)^{nund-1-k} | det N.
#
# For each n and each k in 0..nund-1 we check both:
#   (i)  rank( N|_{s=-(2k+1)} ) <= k+1,
#   (ii) (s+2k+1)^{nund-1-k} divides det N (as a polynomial in s).
#
# Run:  sage check_Lem10_8.sage N

import sys

R = PolynomialRing(QQ, 's'); s = R.gen(); F = R.fraction_field()
def gbin(a, d): return F(0) if d < 0 else prod(a - i for i in range(d)) / factorial(d)

def Nmatrix(nbar, nund):
    b = (s - 3)/2
    return matrix(F, nund, nund,
                  lambda r, m: (2*m+1)*gbin(b, nbar+r-m) + (s+2*m+2)*gbin(b, nbar+r-m-1))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    print("Lemma 10.8 (lem:t3rank) of d05eulershift.tex   [t=3]")
    print("")
    print("  n |  k  | rank N at s=-(2k+1) | <= k+1 ? | (s+2k+1)^(nund-1-k) | div detN ?")
    print("----+-----+---------------------+----------+---------------------+-----------")

    all_ok = True
    for n in range(1, N + 1):
        nbar = (n + 1)//2; nund = n//2
        if nund == 0:
            print("%3d | (nund=0: N is empty; nothing to check)" % n)
            continue
        Nm = Nmatrix(nbar, nund)
        detN = R(Nm.det())
        for k in range(nund):
            rk = matrix(QQ, nund, nund,
                        lambda r, m: QQ(Nm[r, m].subs(s=-(2*k+1)))).rank()
            e = nund - 1 - k
            rank_ok = (rk <= k + 1)
            div_ok = (detN % (s + 2*k + 1)^e == 0) if e >= 1 else True
            if not (rank_ok and div_ok): all_ok = False
            print("%3d | %3d | %19d | %8s | %19d | %s"
                  % (n, k, rk, "yes" if rank_ok else "NO", e, "yes" if div_ok else "NO"))

    print("")
    if all_ok:
        print("Resume: all rank/divisibility checks pass for n = 1..%d. Lemma checked." % N)
    else:
        print("Resume: some check FAILED!")
