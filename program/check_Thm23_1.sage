# check_Thm23_1.sage -- verification of Theorem 23.1 (thm:besselst) of d24Bessel.tex.
#
# The Bessel analogue of the Euler (s,t) family.  With the Bessel cosine
# cosb_s(x) = Gamma(s+1)(2/x)^s J_s(x) = sum (-1)^k x^{2k}/(4^k k! (s+1)_k),
#   f_{s,t}(x) = cosb_s(x) + int_0^x cosb_t(y) dy,
# the coefficient sequence a=(a_k) (a_k = k![x^k] f) is
#   a_{2k}   = (-1)^k (1/2)_k / (s+1)_k,
#   a_{2k+1} = (-1)^k (1/2)_k / (t+1)_k.
# Theorem 23.1 gives the closed form of the dilated Hankel determinant
# HH_n(f) = det(a_{2i+j})_{0<=i,j<n}, with nb=ceil(n/2), nn=floor(n/2):
#   HH_n = (-1)^C(nb,2) * (prod_{i<n}(2i)!)/2^{n(n-1)}
#          * prod_{l=0}^{nb-1} (s+1/2)_l/(s+1)_{n-1+l}
#          * prod_{m=0}^{nn-1} (t+1/2)_m/(t+1)_{n-1+m}
#          * prod_{r=1}^{nb} prod_{c=1}^{nn} (t-s+c-r)/(n-r-c+1).
#
# Both sides lie in QQ(s,t); s,t are kept SYMBOLIC (indeterminates), so each
# LHS(n)=RHS(n) is checked as an identity of rational functions in (s,t).
#
# Run:  sage check_Thm23_1.sage N

import sys

Rst = PolynomialRing(QQ, 's,t'); F = Rst.fraction_field()
s, t = F(Rst.gen(0)), F(Rst.gen(1))

def poch(x, k): return prod(x + j for j in range(k))
H = QQ(1) / 2

def a_moment(midx):
    if midx % 2 == 0:
        k = midx // 2
        return (-1) ** k * poch(H, k) / poch(s + 1, k)
    k = (midx - 1) // 2
    return (-1) ** k * poch(H, k) / poch(t + 1, k)

def LHS(n):
    return matrix(F, n, n, lambda i, j: a_moment(2 * i + j)).det()

def RHS(n):
    nb = (n + 1) // 2; nn = n // 2
    val = (-1) ** binomial(nb, 2) * prod(factorial(2 * i) for i in range(n)) / QQ(2) ** (n * (n - 1))
    val *= prod(poch(s + H, l) / poch(s + 1, n - 1 + l) for l in range(nb))
    val *= prod(poch(t + H, m) / poch(t + 1, n - 1 + m) for m in range(nn))
    val *= prod((t - s + c - r) / (n - r - c + 1)
                for r in range(1, nb + 1) for c in range(1, nn + 1))
    return val

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("Theorem 23.1 (thm:besselst) of d24Bessel.tex   [s,t kept symbolic]")
    print("HH_n(f_{s,t}) = det(a_{2i+j}) = closed product  (identity in QQ(s,t))")
    print("")
    print("  n |  dif = LHS(n) - RHS(n)")
    print("----+-----------------------")

    all_zero = True
    for n in range(1, N + 1):
        d = F(LHS(n) - RHS(n))
        if d != 0:
            all_zero = False
        print("%3d | %s" % (n, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
