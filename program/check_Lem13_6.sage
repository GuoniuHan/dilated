# check_Lem13_6.sage -- verification of Lemma 13.6 (lem:deg) of d07xcosshift.tex.
#
# Degree bound for the shifted determinant.  For f = (1+x)/cos(x)^{s+1} with
#   a_{2k} = E^{(s+1)}_{2k},  a_{2k+1} = (2k+1) E^{(s+1)}_{2k},
# and HH_n^{(1)} = det(a_{2i+j+1})_{0<=i,j<n}, the lemma asserts
#   deg_s HH_n^{(1)} <= C(n,2).
#
# s is kept SYMBOLIC over QQ.  We compute HH_n^{(1)} as a polynomial in s and
# compare its degree with the bound C(n,2).  (In fact equality holds, by
# Theorem 13.1; the lemma only claims "<=".)
#
# Run:  sage check_Lem13_6.sage N
# prints, for n = 1..N, deg_s HH_n^{(1)}, the bound C(n,2) and OK/FAIL,
# then a resume.

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()

def build_a(M):
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2)
    x = Rx.gen()
    g = (1 + x) * (-(s + 1) * cos(x).log()).exp()
    return [factorial(m) * g[m] for m in range(M + 1)]

def HH_shift(n, a):
    return matrix(Rs, n, n, lambda i, j: a[2 * i + j + 1]).det()

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)

    print("Lemma 13.6 (lem:deg) of d07xcosshift.tex   [s kept symbolic]")
    print("deg_s HH_n^(1) <= C(n,2)")
    print("")
    print("  n | deg_s HH_n^(1) | bound C(n,2) | status")
    print("----+----------------+--------------+-------")

    all_ok = True
    for n in range(1, N + 1):
        H = HH_shift(n, a)
        dg = H.degree()
        bd = binomial(n, 2)
        ok = dg <= bd
        if not ok:
            all_ok = False
        print("%3d | %14d | %12d | %s" % (n, dg, bd, "OK" if ok else "FAIL"))

    print("")
    if all_ok:
        print("Resume: all degrees within bound for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some degree exceeds bound. Formula NOT checked!")
