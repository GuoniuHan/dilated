# check_Prop22_1.sage -- verification of Proposition 22.1 (prop:algshift) of d19algesq.tex.
#
# All shifts of the algebraic family f_1(x) = (1+x)(1-x^2)^{-s/2}, with moments
#   a_{2k} = mu_k = (2k)!/k! (s/2)_k,   a_{2k+1} = (2k+1) mu_k.
# The shifted dilated determinant is HH_n^{(r)}(f_1) = det(a_{2i+j+r})_{0<=i,j<n}.
# The proposition claims, for all n >= 1 and r >= 0,
#   HH_n^{(r)}(f_1) = 2^C(n,2) prod_{k=1}^{n-1} k! * prod_{i=0}^{n-1} a_{2i+r}(s),
# equivalently prod_i a_{2i+r} = prod_i mu_{i+rho}  (r even)
#           or prod_i (2i+r) mu_{i+rho}             (r odd),   rho = floor(r/2).
#
# STRONGEST check: s is kept as an INDETERMINATE; LHS and RHS are compared as
# polynomials in s over QQ, for several shifts r.
#
# Run:  sage check_Prop22_1.sage N    (checks n=1..N, r=0..4)

import sys

R = PolynomialRing(QQ, 's')
s = R.gen()

def poch(base, i):
    return prod(base + j for j in range(i))

def mu(k):
    if k < 0:
        return R(0)
    return R(factorial(2 * k) / factorial(k)) * poch(s / 2, k)

def a(m):
    k = m // 2
    return mu(k) if m % 2 == 0 else m * mu(k)      # a_{2k+1}=(2k+1)mu_k, m=2k+1

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    RMAX = 4

    def LHS(n, r):
        return matrix(R, n, n, lambda i, j: a(2 * i + j + r)).det()

    def RHS(n, r):
        return (QQ(2) ** binomial(n, 2) * prod(factorial(k) for k in range(1, n))
                * prod(a(2 * i + r) for i in range(n)))

    def RHS_split(n, r):
        # parity-split product form of prod_i a_{2i+r}
        rho = r // 2
        V = QQ(2) ** binomial(n, 2) * prod(factorial(k) for k in range(1, n))
        if r % 2 == 0:
            pr = prod(mu(i + rho) for i in range(n))
        else:
            pr = prod((2 * i + r) * mu(i + rho) for i in range(n))
        return V * pr

    print("Proposition 22.1 (prop:algshift) of d19algesq.tex   [exact, s indeterminate]")
    print("")
    print("(n,r) | deg | dif = LHS-RHS | dif_split (both should be 0)")
    print("------+-----+---------------+-----------------------------")

    all_zero = True
    for n in range(1, N + 1):
        for r in range(RMAX + 1):
            L = LHS(n, r)
            d = L - RHS(n, r)
            d2 = L - RHS_split(n, r)
            if d != 0 or d2 != 0:
                all_zero = False
            print("(%d,%d) | %3s | %13s | %s" % (n, r, L.degree(), d, d2))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n=1..%d, r=0..%d. Formula checked." % (N, RMAX))
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
