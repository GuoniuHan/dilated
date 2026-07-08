# check_Eq10_8.sage -- verification of identity (10.8) (eq:t3CV) of d05eulershift.tex.
#
# In the rank analysis of N at s=-(2k+1) one writes N_{r,m}=(-1)^u R(r,m),
#   R(r,m) = (2m+1) C(k+u,k) + 2k C(k+u,k+1),   u = X - m,  X = nbar+r,
# and the Chu-Vandermonde identity (10.8) rewrites R as
#   R = sum_{j=0}^{k+1} C(k+X, j) G_j(m),
#       G_j(m) = (2m+1) C(-m, k-j) + 2k C(-m, k+1-j).
#
# Part A checks the identity (10.8) as an integer identity over ranges of (k,X,m).
# Part B checks the link N_{r,m}|_{s=-(2k+1)} = (-1)^u R(r,m) on the actual matrix N.
#
# Run:  sage check_Eq10_8.sage N

import sys

R_ = PolynomialRing(QQ, 's'); s = R_.gen(); F = R_.fraction_field()
def gbin(a, d): return F(0) if d < 0 else prod(a - i for i in range(d)) / factorial(d)

def Rfun(k, X, m):
    u = X - m
    return (2*m+1)*binomial(k+u, k) + 2*k*binomial(k+u, k+1)

def Gj(k, j, m):
    return (2*m+1)*binomial(-m, k-j) + 2*k*binomial(-m, k+1-j)

def CV_rhs(k, X, m):
    return sum(binomial(k+X, j) * Gj(k, j, m) for j in range(k+2))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    all_zero = True

    print("Identity (10.8) (eq:t3CV) of d05eulershift.tex   [t=3]")
    print("")
    print("Part A -- R(k,X,m) = sum_j C(k+X,j) G_j(m)   over 0<=k<=%d, m<=%d, u=X-m in 0..%d"
          % (N, 2*N, 2*N))
    worst = 0
    npairs = 0
    for k in range(N + 1):
        for m in range(2*N + 1):
            for u in range(2*N + 1):
                X = m + u
                d = Rfun(k, X, m) - CV_rhs(k, X, m)
                npairs += 1
                if d != 0: all_zero = False; worst += 1
    print("  checked %d triples (k,X,m); mismatches: %d" % (npairs, worst))

    print("")
    print("Part B -- N_{r,m}|_{s=-(2k+1)} = (-1)^u R(r,m)  on the matrix N (per n,k)")
    print("  n |  k  | max |N - (-1)^u R| over r,m")
    print("----+-----+-----------------------------")
    for n in range(1, N + 1):
        nbar = (n + 1)//2; nund = n//2
        for k in range(nund):
            b = (s - 3)/2
            bad = 0
            for r in range(nund):
                for m in range(nund):
                    u = nbar + r - m
                    Nrm = ((2*m+1)*gbin(b, u) + (s+2*m+2)*gbin(b, u-1)).subs(s=-(2*k+1))
                    val = (-1)^u * Rfun(k, nbar + r, m)
                    if F(Nrm) - val != 0: bad += 1
            if bad: all_zero = False
            print("%3d | %3d | %s" % (n, k, "0" if bad == 0 else "MISMATCH(%d)" % bad))

    print("")
    if all_zero:
        print("Resume: all dif = 0. Identity checked.")
    else:
        print("Resume: some dif != 0. Identity NOT checked!")
