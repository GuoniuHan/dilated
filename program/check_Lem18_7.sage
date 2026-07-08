# check_Lem18_7.sage -- verification of Lemma 18.7 (lem:xsin-split) of d12xsinx.tex.
#
# Catalan splitting of the connection coefficients.  With C_k the Catalan numbers
# and gamma(k) = (2k)!/(4^k k!), for i > m:
#   kappa_{i,m} = alpha_i beta_m C_{i-m-1} C_{i+m},
#     alpha_i = 2(-1)^i/16^i * (2i+1)!(i!)^3 / (gamma(i+1) gamma(2i)),
#     beta_m  = (-1)^{m+1} gamma(m+1)/(m!)^3.
# kappa is taken from its closed form (Lemma 18.6, eq:xsin-kappaform).
#
# Run:  sage check_Lem18_7.sage N   (checks all pairs 0 <= m < i <= N)

import sys

def poch(x, k): return prod(x + j for j in range(k))
def gam(k):     return factorial(2*k) / (4^k * factorial(k))
def Cat(k):     return binomial(2*k, k) / (k + 1)

def kappa_formula(i, m):
    d = i - m
    return (poch(QQ(3)/2 - d, d) * poch(i+m+2, d) * poch(m+1, d)^3
            / (factorial(d) * poch(m + QQ(3)/2, d) * poch(i+m+QQ(1)/2, d)))

def alpha(i): return 2*(-1)^i / 16^i * factorial(2*i+1) * factorial(i)^3 / (gam(i+1)*gam(2*i))
def beta(m):  return (-1)^(m+1) * gam(m+1) / factorial(m)^3

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    print("Lemma 18.7 (lem:xsin-split) of d12xsinx.tex")
    print("checking all pairs 0 <= m < i <= %d" % N)
    print("")
    print("(i,m) |          kappa (closed form)          |  dif = kappa - alpha_i beta_m C_{i-m-1} C_{i+m}")
    print("------+---------------------------------------+-----------------------------------------------")

    all_zero = True
    for i in range(1, N + 1):
        for m in range(i):                       # m < i
            L = kappa_formula(i, m)
            Rhs = alpha(i) * beta(m) * Cat(i - m - 1) * Cat(i + m)
            d = L - Rhs
            if d != 0: all_zero = False
            print("(%d,%d) | %37s | %s" % (i, m, L, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= m < i <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
