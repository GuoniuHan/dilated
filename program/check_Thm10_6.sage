# check_Thm10_6.sage -- verification of Theorem 10.6 (thm:t3) of d05eulershift.tex.
#
# Single shift at t=3.  HH_n^{(1)}|_{t=3} = det(a_{2i+j+1}) factors completely into
# linear forms over Q; all factors have integer roots except one, the carrier
#   Gamma_n(s) = (2nbar+1) s + (4 nbar nund - 2 nbar - 1)        (eq:t3nonclassical)
# whose root is never an integer.  Equivalently (eq:t3collapse),
#   sum_p c_p D_p = A * D_0 * prod_{j=0}^{nund-2}(s+2j+1) * Gamma_n(s),   A>0.
#
# Checks per n:
#   (1) HH_n^{(1)} factors into linear forms over Q  (every irreducible factor deg 1);
#   (2) Gamma_n divides HH_n^{(1)}, its root is non-integer, and it is the UNIQUE
#       factor with non-integer root;
#   (3) collapse eq:t3collapse holds (dif = 0), with A>0.
#
# Run:  sage check_Thm10_6.sage N

import sys

R = PolynomialRing(QQ, 's'); s = R.gen(); F = R.fraction_field()

def gbin(a, d): return F(0) if d < 0 else prod(a - i for i in range(d)) / factorial(d)
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

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    muS = moments(cS, lamS, 2*N + 2)
    muT = moments(cT, lamT, 2*N + 2)
    def a(m): return muS[m//2] if m % 2 == 0 else muT[m//2]

    print("Theorem 10.6 (thm:t3) of d05eulershift.tex   [s symbolic, t=3]")
    print("")
    print("  n | all factors  | carrier Gamma_n(s)   | Gamma_n |  unique   | collapse | A>0")
    print("    | linear/Q ?   | (root non-integer)   | HH ?    | non-int.? | dif=0 ?  |")
    print("----+--------------+----------------------+---------+-----------+----------+----")

    all_ok = True
    for n in range(1, N + 1):
        nbar = (n + 1)//2; nund = n//2
        HH = R(matrix(F, n, n, lambda i, j: a(2*i + j + 1)).det())
        Gam = (2*nbar+1)*s + (4*nbar*nund - 2*nbar - 1)          # carrier
        rc = QQ(-(4*nbar*nund - 2*nbar - 1)) / (2*nbar+1)        # its root

        facs = [f for f, e in HH.factor()]
        linear = all(f.degree() == 1 for f in facs)
        # roots (each linear factor a*s+b -> root -b/a)
        roots = [QQ(-f[0]/f[1]) for f in facs if f.degree() == 1]
        nonint = [r for r in roots if r not in ZZ]
        if nund == 0:
            # HH_n = 1: no carrier exists; carrier claim is vacuous
            divides = True; carrier_root_nonint = True; unique_nonint = True
        else:
            divides = (HH % Gam == 0)
            carrier_root_nonint = (rc not in ZZ)
            unique_nonint = (len(nonint) == 1 and nonint[0] == rc)

        # collapse eq:t3collapse
        b = (s - 3)/2
        Sig = F(0)
        for p in range(nund+1):
            cp = prod(2*m+1 for m in range(p)) * prod(s+2*m+2 for m in range(p, nund))
            cols = [k for k in range(nund+1) if k != p]
            Dp = matrix(F, nund, nund, lambda r, jj: gbin(b, nbar+r-cols[jj])).det()
            Sig += cp*Dp
        if nund == 0:
            collapse_zero = True; Apos = True
        else:
            cols0 = list(range(1, nund+1))
            D0 = matrix(F, nund, nund, lambda r, jj: gbin(b, nbar+r-cols0[jj])).det()
            rhs0 = R(D0 * prod(s+2*j+1 for j in range(nund-1)) * Gam)
            Sig = R(Sig)
            A = Sig.leading_coefficient() / rhs0.leading_coefficient()
            collapse_zero = (Sig - A*rhs0 == 0)
            Apos = (A > 0)

        ok = linear and divides and carrier_root_nonint and unique_nonint and collapse_zero and Apos
        if not ok: all_ok = False
        print("%3d | %12s | %20s | %7s | %9s | %8s | %s"
              % (n, "yes" if linear else "NO", str(Gam) + (" *" if carrier_root_nonint else ""),
                 "yes" if divides else "NO", "yes" if unique_nonint else "NO",
                 "yes" if collapse_zero else "NO", "yes" if Apos else "NO"))

    print("")
    if all_ok:
        print("Resume: all checks pass for n = 1..%d. Theorem checked." % N)
    else:
        print("Resume: some check FAILED!")
