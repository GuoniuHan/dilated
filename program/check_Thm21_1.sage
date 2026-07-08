# check_Thm21_1.sage -- verification of Theorem 21.1 (conj:algsq) of d19algesq.tex.
#
# The squared algebraic family f(x) = (1+x)^2 / (1-x^2)^{s/2}, with moments
#   b_m = m! [x^m] f = m! [x^m] (1+x)^2 exp(-(s/2) log(1-x^2)).
# With mu_k(s) = (2k)!/k! (s/2)_k the moments of (1+x)(1-x^2)^{-s/2}, one has
#   b_{2k}   = mu_k + 2k(2k-1) mu_{k-1},   b_{2k+1} = 2(2k+1) mu_k.
# The theorem asserts, for all n >= 1,
#   HH_n(f) = det(b_{2i+j})_{0<=i,j<n}
#           = 2^C(n,2) prod_{k=1}^{n-1}(2k)! * prod_{i=1}^{n-2}(s/2)_i
#              * [ (2n-1)!! + (s-4) 2^{n-1} (s/2)_{n-1} ] / (s-3).
# The last factor is a monic polynomial in s of degree n-1 (its numerator
# vanishes at s=3), so after clearing (s-3) the identity is a polynomial
# identity in s valid for EVERY s.
#
# The exponent s is kept SYMBOLIC over QQ.  The (s-3) denominator is cleared by
# exact polynomial division (its remainder is checked to be 0).  Two checks:
#   (A) b_m from the EGF matches the mu-formulas b_{2k}, b_{2k+1};
#   (B) det(b_{2i+j}) equals the (denominator-cleared) closed form.
#
# Run:  sage check_Thm21_1.sage N
# prints, for n = 1..N, dif = LHS(n)-RHS(n) (and the cleared-remainder check),
# then a resume.

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()

def poch(a, k):
    return prod(a + j for j in range(k))          # Pochhammer (a)_k

def df(k):
    return prod(range(k, 0, -2))                  # double factorial k!!

def mu(k):
    # mu_k(s) = (2k)!/k! (s/2)_k   (0 for k < 0)
    if k < 0:
        return Rs(0)
    return factorial(2*k) / factorial(k) * poch(s/2, k)

def build_b(M):
    # b_0..b_M from the EGF (1+x)^2 (1-x^2)^{-s/2}
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2)
    x = Rx.gen()
    f = (1 + x)**2 * (-(s/2) * (1 - x**2).log()).exp()
    return [factorial(m) * f[m] for m in range(M + 1)]

def lastfactor(n):
    # [ (2n-1)!! + (s-4) 2^{n-1} (s/2)_{n-1} ] / (s-3), exact division
    num = Rs(df(2*n - 1) + (s - 4) * 2**(n - 1) * poch(s/2, n - 1))
    q, r = num.quo_rem(s - 3)
    return q, r

def RHS(n):
    q, _ = lastfactor(n)
    return (2 ** binomial(n, 2)
            * prod(factorial(2*k) for k in range(1, n))
            * prod(poch(s/2, i) for i in range(1, n - 1))
            * q)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    b = build_b(3 * N)          # max index needed is 3n-3

    print("Theorem 21.1 (conj:algsq) of d19algesq.tex   [s kept symbolic]")
    print("")

    all_ok = True

    print("(A) b_m from the EGF vs the mu-formulas  b_{2k}=mu_k+2k(2k-1)mu_{k-1}, b_{2k+1}=2(2k+1)mu_k")
    print("  k | dif b_{2k} | dif b_{2k+1}")
    print("----+------------+-------------")
    for k in range((3 * N) // 2 + 1):
        d0 = Rs(b[2*k] - (mu(k) + 2*k*(2*k - 1) * mu(k - 1)))
        d1 = Rs(b[2*k + 1] - 2*(2*k + 1) * mu(k)) if 2*k + 1 <= 3*N else Rs(0)
        if d0 != 0 or d1 != 0:
            all_ok = False
        print("%3d | %10s | %s" % (k, d0, d1))

    print("")
    print("(B) HH_n(f) = det(b_{2i+j}) vs the closed form (denominator (s-3) cleared)")
    print("  n | deg last factor | (s-3)-remainder | dif = LHS-RHS")
    print("----+-----------------+-----------------+--------------")
    for n in range(1, N + 1):
        L = Rs(matrix(Rs, n, n, lambda i, j: b[2*i + j]).det())
        q, r = lastfactor(n)
        R = Rs(RHS(n))
        d = L - R
        if d != 0 or r != 0:
            all_ok = False
        print("%3d | %15s | %15s | %s" % (n, q.degree(), r, d))

    print("")
    if all_ok:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
