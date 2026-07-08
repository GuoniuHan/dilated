# check_Lem18_6.sage -- verification of Lemma 18.6 (lem:xsin-kappaclosed) of d12xsinx.tex.
#
# Connection coefficients P_i = sum_m kappa_{i,m} rho_m between the even functional
# S (S-fraction u_n = n^4/(4n^2-1)) and the odd Wilson functional T.  They are the
# unique solution of the connection recurrence (eq:xsin-kapparec)
#   k_{i+1,m} = k_{i,m-1} + (cT_m - cS_i) k_{i,m} + lamT_{m+1} k_{i,m+1} - lamS_i k_{i-1,m}
# with k_{0,0}=1.  Lemma 18.6 (eq:xsin-kappaform):
#   kappa_{i,m} = (3/2-d)_d (i+m+2)_d ((m+1)_d)^3 / [ d! (m+3/2)_d (i+m+1/2)_d ],  d=i-m.
#
# We GENERATE kappa from the recurrence and compare to the closed formula.
# Exact rational arithmetic.
#
# Run:  sage check_Lem18_6.sage N   (checks all pairs 0 <= m <= i <= N)

import sys

def poch(x, k): return prod(x + j for j in range(k))
def uu(n): return QQ(n)^4 / (4*n^2 - 1)

def cS(i):   return uu(2*i) + uu(2*i+1)
def lamS(i): return uu(2*i-1) * uu(2*i)
def cT(m):   return 2*m^2 + 2*m + 1
def lamT(m): return QQ(4*m^6) / ((2*m-1)*(2*m+1)) if m >= 1 else QQ(0)

def kappa_formula(i, m):
    d = i - m
    return (poch(QQ(3)/2 - d, d) * poch(i+m+2, d) * poch(m+1, d)^3
            / (factorial(d) * poch(m + QQ(3)/2, d) * poch(i+m+QQ(1)/2, d)))

def generate_kappa(Imax):
    # rows k[i] as dict m->value; k[i][m]=0 if absent
    K = [dict() for _ in range(Imax + 1)]
    K[0][0] = QQ(1)
    def g(i, m):
        return K[i].get(m, QQ(0)) if 0 <= i <= Imax else QQ(0)
    for i in range(Imax):
        for m in range(0, i + 2):
            val = (g(i, m-1) + (cT(m) - cS(i))*g(i, m)
                   + lamT(m+1)*g(i, m+1) - lamS(i)*g(i-1, m))
            if val != 0: K[i+1][m] = val
    return K

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    K = generate_kappa(N)

    print("Lemma 18.6 (lem:xsin-kappaclosed) of d12xsinx.tex")
    print("checking all pairs 0 <= m <= i <= %d" % N)
    print("")
    print("(i,m) |     kappa (from recurrence)     |  dif = recurrence - formula")
    print("------+---------------------------------+---------------------------")

    all_zero = True
    for i in range(N + 1):
        for m in range(i + 1):
            gen = K[i].get(m, QQ(0))
            d = gen - kappa_formula(i, m)
            if d != 0: all_zero = False
            print("(%d,%d) | %31s | %s" % (i, m, gen, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= m <= i <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
