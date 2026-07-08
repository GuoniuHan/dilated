# check_Thm23_7.sage -- verification of Theorem 23.7 (thm:bessel-det) of d24Bessel.tex.
#
# The connection determinant of the Bessel (s,t) family, evaluated by
# Desnanot--Jacobi condensation.  With the connection coefficients of
# Lemma 23.6 (delta=t-s, d=i-m),
#   kappa_{i,m} = (2i)!/((2m)! 4^d) binom(delta,d) / ((t+2m+1)_d (s+i+m)_d),
# Theorem 23.7 asserts, for integers p>=a>=0, N>=0 and q=p-a,
#   det(kappa_{p+r, a+m})_{0<=r,m<=N-1}
#     = (prod_{r=0}^{N-1} (2p+2r)!/(2a+2r)!) * 4^{-qN}
#       * prod_{r=1}^{q} prod_{c=1}^{N}
#           (delta+c-r) / ( [(q-r)+(N-c)+1] (s+2p-1+c-r) (t+2a+2N-1+r-c) ).
#
# s,t kept SYMBOLIC (identity in QQ(s,t)); checked over the triples (p,a,N)
# with 0<=a<=p<=N and 1<=N<=Nmax (covers every shape of the offset q=p-a).
#
# Run:  sage check_Thm23_7.sage Nmax

import sys

Rst = PolynomialRing(QQ, 's,t'); F = Rst.fraction_field()
s, t = F(Rst.gen(0)), F(Rst.gen(1))
delta = t - s

def poch(x, k): return prod(x + j for j in range(k))
def binomf(x, k): return prod(x - j for j in range(k)) / factorial(k)

def kappa(i, m):
    d = i - m
    if d < 0 or m > i: return F(0)
    return (factorial(2 * i) / factorial(2 * m) / QQ(4) ** d
            * binomf(delta, d) / (poch(t + 2 * m + 1, d) * poch(s + i + m, d)))

def LHS(p, a, N):
    return matrix(F, N, N, lambda r, m: kappa(p + r, a + m)).det()

def RHS(p, a, N):
    q = p - a
    val = prod(factorial(2 * p + 2 * r) / factorial(2 * a + 2 * r) for r in range(N)) * QQ(4) ** (-q * N)
    val *= prod((delta + c - r)
                / (((q - r) + (N - c) + 1) * (s + 2 * p - 1 + c - r) * (t + 2 * a + 2 * N - 1 + r - c))
                for r in range(1, q + 1) for c in range(1, N + 1))
    return val

if __name__ == "__main__":
    Nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 4

    print("Theorem 23.7 (thm:bessel-det) of d24Bessel.tex   [s,t symbolic]")
    print("det(kappa_{p+r,a+m})_{0<=r,m<N} = closed product  (identity in QQ(s,t))")
    print("")
    print(" (p,a,N) | q=p-a |  dif = LHS - RHS")
    print("---------+-------+-----------------")

    all_zero = True
    for N in range(1, Nmax + 1):
        for p in range(0, N + 1):
            for a in range(0, p + 1):
                d = F(LHS(p, a, N) - RHS(p, a, N))
                if d != 0: all_zero = False
                print("(%d,%d,%d) | %5d | %s" % (p, a, N, p - a, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for triples with N=1..%d. Formula checked." % Nmax)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
