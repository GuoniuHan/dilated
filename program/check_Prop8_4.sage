# check_Prop8_4.sage -- verification of Proposition 8.4 (prop:dshift) of d04euler.tex.
#
# Double shift of the Euler (s,t) family.  Let a = (a_n) have
#   a_{2k}   = E^{(s+1)}_{2k}   (moments of the S-fraction u_j = j(j+s)),
#   a_{2k+1} = E^{(t+1)}_{2k}   (moments of the S-fraction v_j = j(j+t)),
# i.e. the coefficient sequence of  1/cos(x)^{s+1} + int_0^x dy/cos(y)^{t+1}.
# With nbar = ceil(n/2), nund = floor(n/2), the double-shifted dilated Hankel
# determinant is
#   HH_n^{(2)} = det( a_{2i+j+2} )_{0<=i,j<=n-1}
#     = (2n-1)!! * prod_{k=0}^{nbar-1}(s+2k+1) * prod_{k=0}^{nund-1}(t+2k+1) * HH_n,
# where HH_n = det( a_{2i+j} ) is the unshifted dilated determinant of the SAME a.
#
# Run:  sage check_Prop8_4.sage N
# prints, for n = 1..N, LHS(n)=HH_n^{(2)}, RHS(n) and dif = LHS(n)-RHS(n),
# then a resume stating whether every dif is 0.

import sys

# ---------- parameters (exact rationals; s=0,t=1 gives the Euler numbers) ----------
# NB: for a non-degenerate test keep (t-s)/2 away from the non-negative integers
# (there the closed form has a zero factor and every determinant vanishes).
s = QQ(1)/2
t = QQ(3)/2
# ---------------------------------------------------------------------------------

def sfrac_moments(bfun, K):
    # moments mu_0..mu_K of the Stieltjes S-fraction with coefficients b_j=bfun(j):
    #   sum_k mu_k x^k = 1/(1 - b_1 x/(1 - b_2 x/(1 - ...))).
    PS = PowerSeriesRing(QQ, 'x', default_prec=K + 2)
    x = PS.gen()
    f = PS(1)
    for j in range(K + 1, 0, -1):
        f = (1 - bfun(j) * x * f)^(-1)
    return [f[k] for k in range(K + 1)]

def build_a(Mmax):
    K = Mmax // 2 + 1
    muS = sfrac_moments(lambda j: j*(j + s), K)   # even part, u_j = j(j+s)
    muT = sfrac_moments(lambda j: j*(j + t), K)   # odd part,  v_j = j(j+t)
    return lambda m: muS[m // 2] if m % 2 == 0 else muT[m // 2]

def dfact(k):
    return prod(range(k, 0, -2))   # (2n-1)!! for odd k; empty product = 1

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    a = build_a(3*N)   # largest index used is 2(n-1)+(n-1)+2 = 3n-1 <= 3N-1

    print("Proposition 8.4 (prop:dshift) of d04euler.tex")
    print("parameters: s = %s, t = %s" % (s, t))
    print("")
    print("  n |            LHS(n)            |            RHS(n)            |  dif")
    print("----+-----------------------------+-----------------------------+-----")

    all_zero = True
    for n in range(1, N + 1):
        nbar = (n + 1) // 2      # ceil(n/2)
        nund = n // 2            # floor(n/2)
        HH   = matrix(QQ, n, n, lambda i, j: a(2*i + j)).det()      # unshifted
        LHS  = matrix(QQ, n, n, lambda i, j: a(2*i + j + 2)).det()  # double shift
        RHS  = (dfact(2*n - 1)
                * prod(s + 2*k + 1 for k in range(nbar))
                * prod(t + 2*k + 1 for k in range(nund))
                * HH)
        d = LHS - RHS
        if d != 0:
            all_zero = False
        print("%3d | %27s | %27s | %s" % (n, LHS, RHS, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
