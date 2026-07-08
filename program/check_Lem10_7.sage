# check_Lem10_7.sage -- verification of Lemma 10.7 (lem:t3atomic) of d05eulershift.tex.
#
# At s = -(2k+1) the even part 1/cos^{s+1} = cos^{2k} is a trig polynomial, so the
# even functional S is atomic on the points y_l = -(2l)^2, l=0..k:
#   (eq:t3atomic)  a_{2j}|_{s=-(2k+1)} = sum_{l=0}^k omega_l (-4 l^2)^j,
#     omega_0 = C(2k,k)/4^k,  omega_l = 2 C(2k,k-l)/4^k  (l>=1).
# Consequently (eq:t3atomicdiv)  (s+2k+1)^{nund-k} | HH_n^{(1)}|_{t=3}, 0<=k<=nund-1.
#
# Part A checks eq:t3atomic (a_{2j} at s=-(2k+1)); Part B checks the divisibility.
#
# Run:  sage check_Lem10_7.sage N

import sys

R = PolynomialRing(QQ, 's'); s = R.gen(); F = R.fraction_field()

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

def omega(k, l):
    return binomial(2*k, k)/4^k if l == 0 else 2*binomial(2*k, k-l)/4^k

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    muS = moments(cS, lamS, 2*N + 2)
    muT = moments(cT, lamT, 2*N + 2)
    def a(m): return muS[m//2] if m % 2 == 0 else muT[m//2]

    all_zero = True

    # ---- Part A: eq:t3atomic (atomic moment formula) ----
    print("Lemma 10.7 (lem:t3atomic) of d05eulershift.tex   [t=3]")
    print("")
    print("Part A -- eq:t3atomic:  a_{2j} at s=-(2k+1) vs atomic sum")
    print("(k,j) |   a_{2j}(-(2k+1))   |   atomic sum   |  dif")
    print("------+---------------------+----------------+-----")
    for k in range(N + 1):
        sk = -(2*k + 1)
        for j in range(2*N + 1):
            lhs = QQ(muS[j].subs(s=sk))
            rhs = sum(omega(k, l) * (-4*l^2)^j for l in range(k + 1))
            d = lhs - rhs
            if d != 0: all_zero = False
            print("(%d,%d) | %19s | %14s | %s" % (k, j, lhs, rhs, d))

    # ---- Part B: eq:t3atomicdiv (divisibility) ----
    print("")
    print("Part B -- eq:t3atomicdiv:  (s+2k+1)^(nund-k) | HH_n^{(1)}")
    print("  n |  k  | exponent nund-k | divides ?")
    print("----+-----+-----------------+----------")
    for n in range(1, N + 1):
        nund = n//2
        HH = R(matrix(F, n, n, lambda i, jj: a(2*i + jj + 1)).det())
        for k in range(nund):
            e = nund - k
            ok = (HH % (s + 2*k + 1)^e == 0)
            if not ok: all_zero = False
            print("%3d | %3d | %15d | %s" % (n, k, e, "yes" if ok else "NO"))

    print("")
    if all_zero:
        print("Resume: all dif = 0 / all divisibilities hold. Lemma checked.")
    else:
        print("Resume: some check FAILED!")
