# check_Prop10_5.sage -- verification of Proposition 10.5 (prop:t3eval) of d05eulershift.tex.
#
# Single shift at t=3 (u_j=j(j+s), v_j=j(j+3)).  Closed evaluation
#   HH_n^{(1)}|_{t=3} = det(a_{2i+j+1})
#     = (-1)^{C(nbar,2)} (prod_{l=0}^{nbar-1} hT_l)
#        (prod_{c=0}^{nund-1} hSp_c/(2c+1)!) (prod_{r=0}^{nund-1} (2nbar+2r)!)
#        * sum_{p=0}^{nund} c_p D_p,
# with hT_l=(2l)!(4)_{2l}, hSp_c=(2c+1)(s+2c+1)(2c)!(s+1)_{2c}, and c_p,D_p the
# Cauchy-Binet data (eq:t3cb).  LHS = det(a_{2i+j+1}) computed directly from the
# moments; RHS assembled from the closed form.  s is symbolic (t=3 fixed), so the
# printed factorisations are exact and dif=0 is a polynomial identity in s.
#
# Run:  sage check_Prop10_5.sage N

import sys

R = PolynomialRing(QQ, 's'); s = R.gen(); F = R.fraction_field()
Py = PolynomialRing(F, 'y'); y = Py.gen()

def poch(x, k): return prod(x + i for i in range(k))
def gbin(a, d): return F(0) if d < 0 else prod(a - i for i in range(d)) / factorial(d)
def uu(j): return j*(j + s)

def cS(i):   return 2*i*(2*i+s) + (2*i+1)*(2*i+1+s)
def lamS(i): return (2*i-1)*(2*i-1+s) * 2*i*(2*i+s)
def cT(m):   return 2*m*(2*m+3) + (2*m+1)*(2*m+4)
def lamT(m): return (2*m-1)*(2*m+2) * 2*m*(2*m+3)

def moments(cf, lf, Kmax):
    sz = Kmax + 1
    J = matrix(F, sz, sz)
    for i in range(sz):
        J[i,i] = cf(i)
        if i+1 < sz: J[i,i+1] = F(1); J[i+1,i] = lf(i+1)
    mu = [F(1)]; Jp = identity_matrix(F, sz)
    for k in range(1, Kmax+1): Jp = Jp*J; mu.append(Jp[0,0])
    return mu

def hT(l):  return factorial(2*l)*poch(4, 2*l)
def hSp(c): return (2*c+1)*(s+2*c+1)*factorial(2*c)*poch(s+1, 2*c)

def cp_Dp_sum(nbar, nund):
    b = (s - 3)/2
    tot = F(0)
    for p in range(nund + 1):
        cp = prod(2*m+1 for m in range(p)) * prod(s+2*m+2 for m in range(p, nund))
        cols = [k for k in range(nund + 1) if k != p]
        Dp = matrix(F, nund, nund,
                    lambda r, jj: gbin(b, nbar + r - cols[jj])).det()
        tot += cp * Dp
    return tot

def fac(poly):
    p = R(poly)
    return "0" if p == 0 else str(factor(p))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    muS = moments(cS, lamS, 2*N + 2)
    muT = moments(cT, lamT, 2*N + 2)
    def a(m): return muS[m//2] if m % 2 == 0 else muT[m//2]

    print("Proposition 10.5 (prop:t3eval) of d05eulershift.tex   [s symbolic, t=3]")
    print("LHS(n) = det(a_{2i+j+1}); RHS(n) = closed evaluation. Shown factored.")
    print("")

    all_zero = True
    for n in range(1, N + 1):
        nbar = (n + 1)//2; nund = n//2
        LHS = R(matrix(F, n, n, lambda i, j: a(2*i + j + 1)).det())
        RHS = F((-1)^binomial(nbar, 2)
                * prod(hT(l) for l in range(nbar))
                * prod(hSp(c)/factorial(2*c+1) for c in range(nund))
                * prod(factorial(2*nbar + 2*r) for r in range(nund))
                * cp_Dp_sum(nbar, nund))
        d = R(LHS) - R(RHS)
        if d != 0: all_zero = False
        print("n=%d:  LHS = %s" % (n, fac(LHS)))
        print("      RHS = %s" % fac(RHS))
        print("      dif = %s" % d)
        print("")

    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
