# check_Lem10_1.sage -- verification of Lemma 10.1 (lem:mixedconn3) of d05eulershift.tex.
#
# Single shift on the line t=3.  u_j = j(j+s), v_j = j(j+3).  Let Q_i be the
# T-orthogonal polynomials (v_j = j(j+3)) and P'_m the monic S'-orthogonal
# polynomials (Christoffel transform of S by y).  Writing Q_i = sum_m kt_{i,m} P'_m,
# Lemma 10.1 claims
#   kt_{i,m} = (2i)!/(2m)!   * C((s-1)/2, i-m)
#            + (s+1)(2i)!/(2m+1)! * C((s-3)/2, i-m-1).
#
# We compute kt_{i,m} INDEPENDENTLY as  S'[Q_i P'_m] / S'[P'_m^2]  (S' orthogonality),
# where S'[y^k] = S[y^{k+1}], and compare to the closed formula.  s is symbolic,
# so each identity is checked as a polynomial identity in s.
#
# Run:  sage check_Lem10_1.sage N   (checks all pairs 0 <= m <= i <= N)

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
def cSp(i):  return uu(2*i+1) + uu(2*i+2)
def lamSp(i):return uu(2*i) * uu(2*i+1)

def monic_ops(cf, lf, Nmax):
    P = [Py(1)]
    if Nmax >= 1: P.append(y - cf(0))
    for i in range(1, Nmax): P.append((y-cf(i))*P[i] - lf(i)*P[i-1])
    return P

def moments(cf, lf, Kmax):
    sz = Kmax + 1
    J = matrix(F, sz, sz)
    for i in range(sz):
        J[i,i] = cf(i)
        if i+1 < sz: J[i,i+1] = F(1); J[i+1,i] = lf(i+1)
    mu = [F(1)]; Jp = identity_matrix(F, sz)
    for k in range(1, Kmax+1): Jp = Jp*J; mu.append(Jp[0,0])
    return mu

def applyf(p, mu):
    return sum(p[k]*mu[k] for k in range(p.degree()+1)) if p != 0 else F(0)

def hSp(c): return (2*c+1)*(s+2*c+1)*factorial(2*c)*poch(s+1, 2*c)

def kt_formula(i, m):
    d = i - m
    return (factorial(2*i)/factorial(2*m) * gbin((s-1)/2, d)
            + (s+1)*factorial(2*i)/factorial(2*m+1) * gbin((s-3)/2, d-1))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    muS  = moments(cS, lamS, 2*N + 3)          # S[y^k]
    muSp = muS[1:]                             # S'[y^k] = S[y^{k+1}]
    Q    = monic_ops(cT,  lamT,  N + 1)        # T-orthogonal
    Pp   = monic_ops(cSp, lamSp, N + 1)        # S'-orthogonal

    def kt_op(i, m):
        return applyf(Q[i]*Pp[m], muSp) / hSp(m)

    print("Lemma 10.1 (lem:mixedconn3) of d05eulershift.tex   [s symbolic, t=3]")
    print("checking all pairs 0 <= m <= i <= %d" % N)
    print("")
    print("(i,m) |               kt (from OPs)               |   dif = OP - formula")
    print("------+-------------------------------------------+---------------------")

    all_zero = True
    for i in range(N + 1):
        for m in range(i + 1):
            L = kt_op(i, m)
            Rf = kt_formula(i, m)
            d = F(L - Rf)
            if d != 0: all_zero = False
            print("(%d,%d) | %41s | %s" % (i, m, L, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= m <= i <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
