# check_Thm17_1.sage -- verification of Theorem 17.1 (thm:deriv) of d11springer.tex.
#
# A derivative of the Springer number family at t = 1.  Let
#   f_s(x) = (cos x + sin x)/(cos x - sin x)^s,   a_k(s) = k! [x^k] f_s,
# equivalently a_k = k! [x^k] (cos x + sin x) exp(s * w),  w = -log(cos x - sin x).
# Then the dilated Hankel determinant HH_n = det(a_{2i+j})_{0<=i,j<n} is
#   HH_n(f_s) = 4^C(n,2) (prod_{k=1}^{n-1} k!)
#                 prod_{j=0}^{n-2} [ (s+2j)(s+2j+1) ]^{n-1-j}.
#
# The exponent s is kept SYMBOLIC over QQ, so each LHS(n)=RHS(n) is a polynomial
# identity in s, valid for all s.
#
# Run:  sage check_Thm17_1.sage N
# prints, for n = 1..N, dif = LHS(n)-RHS(n), then a resume.

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()

def build_a(M):
    # a_0..a_M as polynomials in s, from the EGF (cos x + sin x)(cos x - sin x)^{-s}
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2)
    x = Rx.gen()
    w = -(cos(x) - sin(x)).log()               # zero constant term
    fs = (cos(x) + sin(x)) * (s * w).exp()      # (cos x - sin x)^{-s} = exp(s w)
    return [factorial(k) * fs[k] for k in range(M + 1)]

def LHS(n, a):
    return matrix(Rs, n, n, lambda i, j: a[2 * i + j]).det()

def RHS(n):
    return (4 ** binomial(n, 2)
            * prod(factorial(k) for k in range(1, n))
            * prod(((s + 2*j) * (s + 2*j + 1)) ** (n - 1 - j) for j in range(n - 1)))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)          # max index needed is 3n-3

    print("Theorem 17.1 (thm:deriv) of d11springer.tex   [s kept symbolic]")
    print("")
    print("  n | deg_s HH_n | deg_s RHS |  dif")
    print("----+------------+-----------+-----")

    all_zero = True
    for n in range(1, N + 1):
        L = Rs(LHS(n, a))
        R = Rs(RHS(n))
        d = L - R
        if d != 0:
            all_zero = False
        print("%3d | %10s | %9s | %s" % (n, L.degree(), R.degree(), d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
