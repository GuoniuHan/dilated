# check_Prop13_7.sage -- verification of Proposition 13.7 (prop:const) of d07xcosshift.tex.
#
# The ratio HH_n^{(1)}/HH_n(f) is constant in s.  For f = (1+x)/cos(x)^{s+1} with
#   a_{2k} = E^{(s+1)}_{2k},  a_{2k+1} = (2k+1) E^{(s+1)}_{2k},
#   HH_n^{(1)} = det(a_{2i+j+1}),   HH_n(f) = det(a_{2i+j}),
# the proposition asserts HH_n^{(1)} = c_n^{(1)} prod_{i=1}^{n-1}(s+1)_i for a
# constant c_n^{(1)}; equivalently HH_n^{(1)}/HH_n(f) is independent of s.
#
# s is kept SYMBOLIC over QQ.  We form the ratio HH_n^{(1)}/HH_n(f) in QQ(s) and
# verify it lies in QQ (no s-dependence); the constant is printed and compared
# with ((n-1)!!)^2 (its value, by Theorem 13.1).
#
# Run:  sage check_Prop13_7.sage N
# prints, for n = 1..N, the constant ratio, whether it is s-free, and
# dif = ratio - ((n-1)!!)^2, then a resume.

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()
Fs = Rs.fraction_field()

def df(k):
    return prod(range(k, 0, -2))

def build_a(M):
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2)
    x = Rx.gen()
    g = (1 + x) * (-(s + 1) * cos(x).log()).exp()
    return [factorial(m) * g[m] for m in range(M + 1)]

def HH_shift(n, a):
    return matrix(Fs, n, n, lambda i, j: Fs(a[2 * i + j + 1])).det()

def HH(n, a):
    return matrix(Fs, n, n, lambda i, j: Fs(a[2 * i + j])).det()

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)

    print("Proposition 13.7 (prop:const) of d07xcosshift.tex   [s kept symbolic]")
    print("HH_n^(1)/HH_n(f) is constant in s  (= ((n-1)!!)^2)")
    print("")
    print("  n |   ratio   | s-free? | ((n-1)!!)^2 |  dif")
    print("----+-----------+---------+-------------+-----")

    all_ok = True
    for n in range(1, N + 1):
        ratio = HH_shift(n, a) / HH(n, a)
        # s-free iff numerator and denominator are both degree-0 polynomials
        sfree = (ratio.numerator().degree() == 0
                 and ratio.denominator().degree() == 0)
        val = QQ(ratio) if sfree else None
        target = df(n - 1) ** 2
        d = (val - target) if sfree else "n/a"
        if not sfree or d != 0:
            all_ok = False
        print("%3d | %9s | %7s | %11s | %s"
              % (n, val if sfree else "NOT CONST", "yes" if sfree else "NO", target, d))

    print("")
    if all_ok:
        print("Resume: ratio is a constant equal to ((n-1)!!)^2 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: ratio not constant (or != ((n-1)!!)^2). Formula NOT checked!")
