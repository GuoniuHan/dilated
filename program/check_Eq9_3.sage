# check_Eq9_3.sage -- verification of equation (9.3) (eq:shiftconst) of d05eulershift.tex.
#
# The scalar K_n of Proposition 9.1 (single shift, line t=1) is fixed by
# the secant/tangent specialisation s=0, giving eq (9.3):
#   K_n = (2n-1)!! * 2^{-binom(n,2)} * prod_{k=1}^{n-1} k! (2k)!
#              / [ prod_{j=1}^{nn} (1/2)_j^2 * prod_{j=1}^{nn-1} j! * prod_{j=1}^{nb-1} (1/2)_j ],
# with nb = ceil(n/2), nn = floor(n/2), and (1/2)_j = (2j-1)!!/2^j = (2j)!/(4^j j!).
# The first values are stated to be
#   K_n = 1, 12, 1440, 9676800, 7023034368000, 149513097738977280000, ...
#
# This checker verifies (9.3) in THREE independent ways for each n:
#   (i)  the closed formula equals the tabulated first values;
#   (ii) the closed formula equals the s=0 specialisation of the determinant
#        divided by the Pochhammer products, i.e. it really is the scalar
#        K_n of eq (9.2), computed directly from HH_n^{(1)}|_{s=0};
#   (iii) both forms of the denominator (2j-1)!!/2^j and (2j)!/(4^j j!) agree.
# All differences must be 0 (exact rational arithmetic).
#
# Run:  sage check_Eq9_3.sage N

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()

def poch(x, k):
    return prod(x + j for j in range(k))

def build_a(M):
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2)
    x = Rx.gen()
    logcos = cos(x).log()
    Ssec = (-(s + 1) * logcos).exp()
    Tsec = (-2 * logcos).exp()
    a = []
    for m in range(M + 1):
        if m % 2 == 0:
            a.append(factorial(m) * Ssec[m])
        else:
            a.append(factorial(m - 1) * Tsec[m - 1])
    return a

def K_formula(n):
    # eq (9.3) as written, with (1/2)_j via poch().
    nb = (n + 1) // 2
    nn = n // 2
    dfact = prod(2 * k - 1 for k in range(1, n + 1))          # (2n-1)!!
    num = dfact * QQ(2) ** (-binomial(n, 2)) \
          * prod(factorial(k) * factorial(2 * k) for k in range(1, n))
    den = prod(poch(QQ(1) / 2, j) ** 2 for j in range(1, nn + 1)) \
          * prod(factorial(j) for j in range(1, nn)) \
          * prod(poch(QQ(1) / 2, j) for j in range(1, nb))
    return num / den

def K_from_det(n, a):
    # eq (9.2) at s=0:  HH_n^{(1)}|_{s=0} / (Pochhammer products)|_{s=0}.
    nb = (n + 1) // 2
    nn = n // 2
    det0 = matrix(QQ, n, n,
                  lambda i, j: a[2 * i + j + 1].subs(s=0)).det()
    poch0 = prod(poch(QQ(1) / 2, j) ** 2 for j in range(1, nn + 1)) \
          * prod(poch(QQ(1), j) for j in range(1, nn)) \
          * prod(poch(QQ(1) / 2, j) for j in range(1, nb))
    return det0 / poch0

def halfpoch_alt(j):
    # the two claimed closed forms of (1/2)_j
    return (prod(2 * k - 1 for k in range(1, j + 1)) / QQ(2) ** j,
            factorial(2 * j) / (QQ(4) ** j * factorial(j)))

FIRST = [1, 12, 1440, 9676800, 7023034368000, 149513097738977280000]

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)

    print("Equation (9.3) (eq:shiftconst) of d05eulershift.tex   [t=1]")
    print("K_n from formula vs first values vs det|_{s=0}/poch|_{s=0}.")
    print("")
    print("  n |          K_n (formula)           | dif(vals) | dif(det) | dif(1/2poch)")
    print("----+----------------------------------+-----------+----------+-------------")

    all_zero = True
    for n in range(1, N + 1):
        Kf = K_formula(n)
        d_vals = Kf - FIRST[n - 1] if n - 1 < len(FIRST) else QQ(0)
        d_det = Kf - K_from_det(n, a)
        # cross-check both closed forms of (1/2)_j for all j up to n
        d_half = max((abs(p1 - p2) for j in range(1, n + 1)
                      for (p1, p2) in [halfpoch_alt(j)]), default=QQ(0))
        if d_vals != 0 or d_det != 0 or d_half != 0:
            all_zero = False
        print("%3d | %32s | %9s | %8s | %s" % (n, Kf, d_vals, d_det, d_half))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
