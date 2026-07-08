# check_Prop7_1.sage -- verification of Proposition 7.1 (prop:family) of d04euler.tex.
#
# The Euler (secant--tangent) (s,t) family.  Let
#   f(x) = 1/cos(x)^{s+1} + int_0^x dy/cos(y)^{t+1},
# so that the coefficient sequence a = (a_n) is
#   a_{2k}   = E_{2k}^{(s+1)} = (2k)! [x^{2k}] cos(x)^{-(s+1)}   (even part),
#   a_{2k+1} = E_{2k}^{(t+1)} = (2k)! [x^{2k}] cos(x)^{-(t+1)}   (odd part).
# Then the dilated Hankel determinant HH_n(f) = det(a_{2i+j})_{0<=i,j<=n-1} is
#   (-1)^C(nb,2) prod_{i=0}^{n-1}(2i)!
#     * prod_{l=0}^{nb-1}(s+1)_{2l} * prod_{m=0}^{nn-1}(t+1)_{2m}
#     * prod_{r=1}^{nb} prod_{c=1}^{nn} ((t-s)/2 + c - r)/(n - r - c + 1),
# where nb = ceil(n/2), nn = floor(n/2).
#
# Both s and t are kept SYMBOLIC (indeterminates over QQ), so each LHS(n)=RHS(n)
# is checked as a polynomial identity in (s,t), valid for ALL parameter values.
#
# Run:  sage check_Prop7_1.sage N
# prints, for n = 1..N, LHS(n), RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

Rst = PolynomialRing(QQ, 's, t')
s, t = Rst.gens()

def poch(x, k):
    # Pochhammer (x)_k = x(x+1)...(x+k-1)
    return prod(x + j for j in range(k))

def build_a(M):
    # a_0,...,a_M as polynomials in (s,t), from the exponential generating
    # functions cos(x)^{-(s+1)} (even part) and cos(x)^{-(t+1)} (odd part).
    Rx = PowerSeriesRing(Rst, 'x', default_prec=M + 2)
    x = Rx.gen()
    logcos = cos(x).log()                 # -x^2/2 - x^4/12 - ...
    Ssec = (-(s + 1) * logcos).exp()      # sum E_{2k}^{(s+1)} x^{2k}/(2k)!
    Tsec = (-(t + 1) * logcos).exp()      # sum E_{2k}^{(t+1)} x^{2k}/(2k)!
    a = []
    for m in range(M + 1):
        if m % 2 == 0:
            a.append(factorial(m) * Ssec[m])
        else:
            a.append(factorial(m - 1) * Tsec[m - 1])
    return a

def LHS(n, a):
    return matrix(Rst, n, n, lambda i, j: a[2 * i + j]).det()

def RHS(n):
    nb = (n + 1) // 2          # ceil(n/2)
    nn = n // 2                # floor(n/2)
    sign = (-1) ** binomial(nb, 2)
    fac  = prod(factorial(2 * i) for i in range(n))
    sp   = prod(poch(s + 1, 2 * l) for l in range(nb))
    tp   = prod(poch(t + 1, 2 * m) for m in range(nn))
    dp   = prod(((t - s) / 2 + c - r) / (n - r - c + 1)
                for r in range(1, nb + 1) for c in range(1, nn + 1))
    return sign * fac * sp * tp * dp

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)         # max index needed is 2(n-1)+(n-1) = 3n-3

    print("Proposition 7.1 (prop:family) of d04euler.tex   [s, t kept symbolic]")
    print("")
    print("  n |                     LHS(n)                     |               dif")
    print("----+------------------------------------------------+------------------")

    all_zero = True
    for n in range(1, N + 1):
        L = LHS(n, a)
        R = RHS(n)
        d = Rst(L - R)
        if d != 0:
            all_zero = False
        print("%3d | %46s | %s" % (n, L, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
