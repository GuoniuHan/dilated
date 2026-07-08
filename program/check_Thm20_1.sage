# check_Thm20_1.sage -- verification of Theorem 20.1 (conj:alg) of d18alge.tex.
#
# The algebraic family  f(x) = (1+x)/(1-x^2)^{s/2}.  Its coefficient sequence is
#   a_{2k}   = (2k)!/k! * (s/2)_k = 4^k (1/2)_k (s/2)_k,
#   a_{2k+1} = (2k+1) a_{2k}.
# The theorem claims, for all n >= 1,
#   HH_n(f) = ( prod_{k=1}^{n-1} (2k)! ) * prod_{j=1}^{n-1} (s+2j-2)^{n-j}
#           = 2^C(n,2) * prod_{k=1}^{n-1} (2k)! * prod_{i=1}^{n-1} (s/2)_i.
#
# STRONGEST check: s is kept as an INDETERMINATE, so LHS = det(a_{2i+j}) and BOTH
# right-hand forms are compared as polynomials in s over QQ (dif is the zero poly).
#
# Run:  sage check_Thm20_1.sage N
# prints, for n=1..N, dif1 (form 1), dif2 (form 2) and LHS factored, then a resume.

import sys

R = PolynomialRing(QQ, 's')
s = R.gen()

def poch(base, i):
    # rising factorial (base)_i = base (base+1) ... (base+i-1)
    return prod(base + j for j in range(i))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    def a(m):
        k = m // 2
        even = R(factorial(2 * k) / factorial(k)) * poch(s / 2, k)   # a_{2k}
        return even if m % 2 == 0 else m * even                       # a_{2k+1}=(2k+1)a_{2k}, m=2k+1

    def LHS(n):
        return matrix(R, n, n, lambda i, j: a(2 * i + j)).det()

    def cn(n):
        return prod(factorial(2 * k) for k in range(1, n))

    def RHS1(n):
        return R(cn(n)) * prod((s + 2 * j - 2) ** (n - j) for j in range(1, n))

    def RHS2(n):
        return QQ(2) ** binomial(n, 2) * R(cn(n)) * prod(poch(s / 2, i) for i in range(1, n))

    print("Theorem 20.1 (conj:alg) of d18alge.tex   [exact, s indeterminate over QQ]")
    print("a = %s ..." % [a(m) for m in range(6)])
    print("")
    print("  n | deg | dif1=LHS-RHS1 | dif2=LHS-RHS2 (both 0);  LHS factored")
    print("----+-----+---------------+----------------------------------------")

    all_zero = True
    for n in range(1, N + 1):
        L = LHS(n)
        d1 = L - RHS1(n)
        d2 = L - RHS2(n)
        if d1 != 0 or d2 != 0:
            all_zero = False
        fac = L.factor() if L != 0 else "0"
        print("%3d | %3s | %13s | %s ;  LHS=%s" % (n, L.degree(), d1, d2, fac))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
