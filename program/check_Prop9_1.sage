# check_Prop9_1.sage -- verification of Proposition 9.1 (prop:shift) of d05eulershift.tex.
#
# The single shift of the Euler (s,t) family on the line t=1.  With
#   u_j = j(j+s)   (even part = 1/cos(x)^{s+1}),
#   v_j = j(j+1)   (odd part  = the tangent fraction, t=1),
# the sequence a=(a_n) has
#   a_{2k}   = E_{2k}^{(s+1)} = (2k)! [x^{2k}] cos(x)^{-(s+1)}   (even part),
#   a_{2k+1} = E_{2k}^{(2)}   = (2k)! [x^{2k}] cos(x)^{-2}       (odd part).
# The SINGLE-shifted dilated determinant is
#   HH_n^{(1)} = det(a_{2i+j+1})_{0<=i,j<=n-1}.
# Proposition 9.1, eq (9.2), asserts the closed form
#   HH_n^{(1)} = Ktilde_n * prod_{j=1}^{nn} ((s+1)/2)_j^2
#                        * prod_{j=1}^{nn-1} (s/2+1)_j
#                        * prod_{j=1}^{nb-1} ((1-s)/2)_j,
# with nb = ceil(n/2), nn = floor(n/2), and Ktilde_n a positive scalar
# independent of s given by eq (9.4) (eq:shiftconst); it is substituted here.
#
# s is kept SYMBOLIC (indeterminate over QQ), so each LHS(n)=RHS(n) is checked
# as a polynomial identity in s, valid for ALL s (on the line t=1).
#
# Run:  sage check_Prop9_1.sage N
# prints, for n = 1..N, LHS(n), the difference dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()

def poch(x, k):
    # Pochhammer (x)_k = x(x+1)...(x+k-1)
    return prod(x + j for j in range(k))

def build_a(M):
    # a_0,...,a_M as polynomials in s, from the EGFs cos(x)^{-(s+1)} (even part)
    # and cos(x)^{-2} (odd part, t=1).
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2)
    x = Rx.gen()
    logcos = cos(x).log()                 # -x^2/2 - x^4/12 - ...
    Ssec = (-(s + 1) * logcos).exp()      # sum E_{2k}^{(s+1)} x^{2k}/(2k)!
    Tsec = (-2 * logcos).exp()            # sum E_{2k}^{(2)}   x^{2k}/(2k)!
    a = []
    for m in range(M + 1):
        if m % 2 == 0:
            a.append(factorial(m) * Ssec[m])
        else:
            a.append(factorial(m - 1) * Tsec[m - 1])
    return a

def LHS(n, a):
    # single shift: column index advanced by one -> a_{2i+j+1}
    return matrix(Rs, n, n, lambda i, j: a[2 * i + j + 1]).det()

def Ktilde(n):
    # eq (9.4), eq:shiftconst: the positive scalar independent of s.
    nb = (n + 1) // 2
    nn = n // 2
    dfact = prod(2 * k - 1 for k in range(1, n + 1))          # (2n-1)!!
    num = dfact * QQ(2) ** (-binomial(n, 2)) \
          * prod(factorial(k) * factorial(2 * k) for k in range(1, n))
    den = prod(poch(QQ(1) / 2, j) ** 2 for j in range(1, nn + 1)) \
          * prod(factorial(j) for j in range(1, nn)) \
          * prod(poch(QQ(1) / 2, j) for j in range(1, nb))
    return num / den

def RHS(n):
    nb = (n + 1) // 2
    nn = n // 2
    pochpart = prod(poch((s + 1) / 2, j) ** 2 for j in range(1, nn + 1)) \
             * prod(poch(s / 2 + 1, j) for j in range(1, nn)) \
             * prod(poch((1 - s) / 2, j) for j in range(1, nb))
    return Ktilde(n) * pochpart

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)                    # max index needed is 2(n-1)+(n-1)+1 = 3n-2

    print("Proposition 9.1 (prop:shift) of d05eulershift.tex   [s kept symbolic, t=1]")
    print("HH_n^{(1)} = det(a_{2i+j+1}) = Ktilde_n * (Pochhammer products), eq (9.2).")
    print("")
    print("  n |                     LHS(n) = HH_n^{(1)}                      |     dif")
    print("----+--------------------------------------------------------------+--------")

    all_zero = True
    for n in range(1, N + 1):
        L = LHS(n, a)
        R = RHS(n)
        d = Rs(L - R)
        if d != 0:
            all_zero = False
        print("%3d | %60s | %s" % (n, L, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
