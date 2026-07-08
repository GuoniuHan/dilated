# check_Thm16_1.sage -- verification of Theorem 16.1 (thm:springer) of d11springer.tex.
#
# Springer number family.  For every integer r >= 1,
#   HH_n( 1/(cos x - t sin x)^r ) = ( t(t^2+1) )^C(n,2) * HH_n( 1/(1-x)^r ).
# Here HH_n(f) = det(a_{2i+j}); the moments of the left side are
#   a_n(t) = n! [x^n] (cos x - t sin x)^{-r},
# and the right side uses b_n = (r)_n (Pochhammer), the coefficients of 1/(1-x)^r.
#
# t is kept SYMBOLIC, so LHS(n) = RHS(n) is a polynomial identity in t; the LHS is
# printed factored, exhibiting the ( t(t^2+1) )^C(n,2) divisor.  (r=t=1 is the
# Springer case HH_n = 4^C(n,2) prod k!(2k)!.)
#
# Run:  sage check_Thm16_1.sage N

import sys

# ---------- parameter (exponent, integer >= 1) ----------
r = 1
# --------------------------------------------------------

R = PolynomialRing(QQ, 't'); t = R.gen(); F = R.fraction_field()

def springer_moments(Mmax):
    prec = Mmax + 1
    PS = PowerSeriesRing(F, 'x', default_prec=prec + 1)
    x = PS.gen()
    cosx = sum((-1)^k * x^(2*k) / factorial(2*k) for k in range(prec//2 + 1))
    sinx = sum((-1)^k * x^(2*k+1) / factorial(2*k+1) for k in range(prec//2 + 1))
    h = (cosx - t*sinx).add_bigoh(prec + 1)
    f = h.inverse()^r                      # (cos x - t sin x)^{-r}
    return [factorial(m) * f[m] for m in range(Mmax + 1)]   # a_m(t)

def poch(x, k): return prod(x + i for i in range(k))

def fac(poly):
    p = R(poly)
    return "0" if p == 0 else str(factor(p))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = springer_moments(3*N)                     # a_m(t), m up to 3n-3
    def b(m): return poch(r, m)                   # coeffs of 1/(1-x)^r

    print("Theorem 16.1 (thm:springer) of d11springer.tex   [t symbolic]")
    print("parameter: r = %d" % r)
    print("LHS(n) = HH_n(1/(cos x - t sin x)^r) shown factored; RHS = (t(t^2+1))^C(n,2)*HH_n(1/(1-x)^r).")
    print("")

    all_zero = True
    for n in range(1, N + 1):
        LHS = R(matrix(F, n, n, lambda i, j: a[2*i + j]).det())
        c   = matrix(QQ, n, n, lambda i, j: b(2*i + j)).det()   # HH_n(1/(1-x)^r)
        RHS = R((t*(t^2 + 1))^binomial(n, 2) * c)
        d = LHS - RHS
        if d != 0: all_zero = False
        print("n=%d:  HH_n(1/(1-x)^r) = %s" % (n, c))
        print("      LHS = %s" % fac(LHS))
        print("      dif = %s" % d)
        print("")

    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
