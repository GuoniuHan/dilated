# check_Thm18_8.sage -- verification of Theorem 18.8 (thm:xsin-catalan) of d12xsinx.tex.
#
# The Catalan determinant
#   F_N(a,b) = det( C_{a+r-m} C_{b+r+m} )_{0<=r,m<=N-1}
# (Hadamard product of a Toeplitz and a Hankel matrix of Catalan numbers) equals,
# for integers a >= N-1 and b >= a+1, with gamma(k)=(2k)!/(4^k k!) and the
# staircase products Lambda_N(t) = prod_{j=1}^{N-1} (t+j)_j :
#   F_N(a,b) = (-4)^C(N,2) 4^{N(a+b)} prod_{k=1}^{N-1}(2k+1)!
#     * [ prod_{k=0}^{N-1} gamma(a-k) gamma(b+k) ]
#       / [ prod_{k=2}^{N+1}(a+k-1)! * prod_{k=N+1}^{2N}(b+k-1)! ]
#     * Lambda_N(b-a) Lambda_N(a+b+3/2).
#
# For each N we test several in-range (a,b); F is the direct determinant, G the
# closed form.  (The hypothesis b>=a+1 is essential -- see the out-of-range note.)
#
# Run:  sage check_Thm18_8.sage N

import sys

def poch(x, k): return prod(x + j for j in range(k))
def gam(k):     return factorial(2*k) / (4^k * factorial(k))
def Cat(k):     return binomial(2*k, k) / (k + 1)
def Lam(N, t):  return prod(poch(t + j, j) for j in range(1, N))       # Lambda_N(t)

def F(N, a, b):
    if N == 0: return QQ(1)
    return matrix(QQ, N, N, lambda r, m: Cat(a + r - m) * Cat(b + r + m)).det()

def G(N, a, b):
    if N == 0: return QQ(1)
    num = prod(gam(a - k) * gam(b + k) for k in range(N))
    den = prod(factorial(a + k - 1) for k in range(2, N + 2)) * \
          prod(factorial(b + k - 1) for k in range(N + 1, 2*N + 1))
    return ((-4)^binomial(N, 2) * 4^(N*(a + b)) * prod(factorial(2*k+1) for k in range(1, N))
            * num / den * Lam(N, b - a) * Lam(N, a + b + QQ(3)/2))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    print("Theorem 18.8 (thm:xsin-catalan) of d12xsinx.tex")
    print("F_N(a,b)=det(C_{a+r-m} C_{b+r+m}); checked for a>=N-1, b>=a+1")
    print("")
    print(" N  (a, b) |             F_N(a,b)            |  dif = F - G")
    print("-----------+---------------------------------+-------------")

    all_zero = True
    for n in range(0, N + 1):
        for a in range(n - 1 if n >= 1 else 0, n + 2):
            if a < 0: continue
            for b in (a + 1, a + 2):
                Fv = F(n, a, b)
                d = Fv - G(n, a, b)
                if d != 0: all_zero = False
                print("%2d (%2d,%2d) | %31s | %s" % (n, a, b, Fv, d))

    print("")
    # Out-of-range observation (b <= a): report agreement, and the (2,2,2) value.
    ndiff = sum(1 for n in range(2, N + 1) for a in range(n - 1, n + 3)
                for b in range(max(0, a - 2), a + 1) if F(n, a, b) != G(n, a, b))
    print("Out-of-range observation (b <= a): F != G at %d of the scanned points." % ndiff)
    print("  In particular F_2(2,2)=G_2(2,2)=%s (the paper's Remark states the RHS is -13/24" % F(2, 2, 2))
    print("  there; the formula as written gives -13, matching F -- a likely erratum).")
    print("")
    if all_zero:
        print("Resume: all dif = 0 (in range a>=N-1, b>=a+1). Formula checked.")
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
