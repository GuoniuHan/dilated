# check_Lem17_3.sage -- verification of Lemma 17.3 (lem:gausskernel) of d11springer.tex.
#
# Gaussian kernel determinant.  The lemma asserts, as sigma -> 0,
#   det( e^{sigma (2i+j)^2} )_{0<=i,j<n}
#      = 4^C(n,2) (prod_{k=1}^{n-1} k!) sigma^C(n,2) + O(sigma^{C(n,2)+1}).
#
# We compute the determinant as a power series in sigma (exact over QQ) and check
#   (a) its valuation equals C(n,2)  (all lower-order terms vanish), and
#   (b) its coefficient of sigma^C(n,2) equals 4^C(n,2) prod_{k=1}^{n-1} k!.
#
# Run:  sage check_Lem17_3.sage N
# prints, for n = 1..N, the valuation, the leading coefficient and its target,
# with dif = coeff - target, then a resume.

import sys

def build_ring(N):
    return PowerSeriesRing(QQ, 'sig', default_prec=binomial(N, 2) + 3)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    Rsig = build_ring(N)
    sig = Rsig.gen()

    print("Lemma 17.3 (lem:gausskernel) of d11springer.tex   [exact power series over QQ]")
    print("det(e^{sig (2i+j)^2}) = 4^C(n,2) (prod_{k<n} k!) sig^C(n,2) + O(sig^{C(n,2)+1})")
    print("")
    print("  n | val | C(n,2) | leading coeff | 4^C(n,2) prod k! |  dif")
    print("----+-----+--------+---------------+------------------+-----")

    all_ok = True
    for n in range(1, N + 1):
        D = matrix(Rsig, n, n, lambda i, j: (sig * (2*i + j)**2).exp()).det()
        c2 = binomial(n, 2)
        val = D.valuation()
        coeff = D[c2]
        target = 4 ** c2 * prod(factorial(k) for k in range(1, n))
        d = coeff - target
        ok = (val == c2 and d == 0)
        if not ok:
            all_ok = False
        print("%3d | %3d | %6d | %13s | %16s | %s" % (n, val, c2, coeff, target, d))

    print("")
    if all_ok:
        print("Resume: valuation = C(n,2) and leading coeff matches for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some check fails. Formula NOT checked!")
