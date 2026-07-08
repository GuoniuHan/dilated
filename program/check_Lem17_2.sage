# check_Lem17_2.sage -- verification of Lemma 17.2 (lem:Apoly) of d11springer.tex.
#
# Polynomial structure of the moments a_k(s) = k! [x^k] (cos x + sin x)(cos x - sin x)^{-s}.
# The lemma asserts:
#   (i)  a_k(s) is a MONIC polynomial in s of degree k;
#   (ii) its coefficient  pi_r(k) := [s^{k-r}] a_k(s)  is a polynomial in k of
#        degree <= 2r with leading coefficient [k^{2r}] pi_r = 1/r!.
#
# Everything is exact over QQ.  For (ii), pi_r(k) is sampled at k = r,...,r+2r+2
# and Lagrange-interpolated in k; the interpolant's degree (<= 2r ?) and its
# coefficient of k^{2r} (= 1/r! ?) are then checked.
#
# Run:  sage check_Lem17_2.sage N
# prints two tables (monicity; pi_r structure), then a resume.

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()
Rk = PolynomialRing(QQ, 'k')

def build_a(M):
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2)
    x = Rx.gen()
    w = -(cos(x) - sin(x)).log()
    fs = (cos(x) + sin(x)) * (s * w).exp()
    return [factorial(k) * fs[k] for k in range(M + 1)]

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    # need a_k for k up to r + 2r + 2 = 3r+2 with r up to N
    a = build_a(3 * N + 2)

    print("Lemma 17.2 (lem:Apoly) of d11springer.tex   [exact over QQ]")
    print("")

    all_ok = True

    print("(i) a_k(s) is monic of degree k")
    print("  k | deg_s a_k | leading coeff | status")
    print("----+-----------+---------------+-------")
    for k in range(N + 1):
        dg = a[k].degree()
        lc = a[k].leading_coefficient()
        ok = (dg == k and lc == 1)
        if not ok:
            all_ok = False
        print("%3d | %9d | %13s | %s" % (k, dg, lc, "OK" if ok else "FAIL"))

    print("")
    print("(ii) pi_r(k) = [s^{k-r}] a_k(s) is poly in k, deg <= 2r, [k^{2r}] = 1/r!")
    print("  r | deg_k pi_r | bound 2r | [k^{2r}]pi_r |   1/r!   | status")
    print("----+------------+----------+--------------+----------+-------")
    for r in range(N + 1):
        # sample pi_r(k) = [s^{k-r}] a_k for enough k >= r to over-determine deg 2r
        pts = [(k, a[k][k - r]) for k in range(r, r + 2 * r + 3)]
        poly = Rk.lagrange_polynomial(pts)
        dg = poly.degree()
        lead = poly[2 * r]
        target = QQ(1) / factorial(r)
        ok = (dg <= 2 * r and lead == target)
        if not ok:
            all_ok = False
        print("%3d | %10d | %8d | %12s | %8s | %s"
              % (r, dg, 2 * r, lead, target, "OK" if ok else "FAIL"))

    print("")
    if all_ok:
        print("Resume: all checks pass for indices <= %d. Formula checked." % N)
    else:
        print("Resume: some check fails. Formula NOT checked!")
