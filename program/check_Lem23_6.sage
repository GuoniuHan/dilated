# check_Lem23_6.sage -- verification of Lemma 23.6 (lem:bessel-kappa) of d24Bessel.tex.
#
# Connection coefficients of the Bessel (s,t) family.  Let P_i, Q_m be the
# monic orthogonal polynomials of the even functional S (S[y^k]=a_{2k}) and the
# odd functional T (T[y^k]=a_{2k+1}), and P_i = sum_m kappa_{i,m} Q_m the
# connection expansion.  With delta = t-s and d = i-m, Lemma 23.6 asserts
#   kappa_{i,m} = (2i)!/((2m)! 4^d) * binom(delta,d) / ((t+2m+1)_d (s+i+m)_d).
#
# P_i, Q_m are built here from the three-term recurrences whose data come from
# the S-fraction coefficients of Lemma 23.5 (eq:bessel-u),
#   u_j = -j(j+2s-1)/(4(s+j-1)(s+j)),   v_j the same with t,
# via c^S_0=u_1, c^S_i=u_{2i}+u_{2i+1}, lam^S_i=u_{2i-1}u_{2i} (and likewise T).
# LHS(i,m) is read off the expansion of P_i in the Q-basis; RHS(i,m) is the
# closed form.  s,t kept SYMBOLIC (identity in QQ(s,t)); indexed by pairs
# (i,m) with 0<=m<=i<=N.  A built-in orthogonality sanity check confirms
# S[P_i y^j]=0 for j<i.
#
# Run:  sage check_Lem23_6.sage N

import sys

Rst = PolynomialRing(QQ, 's,t'); F = Rst.fraction_field()
s, t = F(Rst.gen(0)), F(Rst.gen(1))
Ry = PolynomialRing(F, 'y'); y = Ry.gen()
delta = t - s
H = QQ(1) / 2

def poch(x, k): return prod(x + j for j in range(k))
def binomf(x, k): return prod(x - j for j in range(k)) / factorial(k)

def a_moment(midx):
    if midx % 2 == 0:
        k = midx // 2
        return (-1) ** k * poch(H, k) / poch(s + 1, k)
    k = (midx - 1) // 2
    return (-1) ** k * poch(H, k) / poch(t + 1, k)

def u(j): return -j * (j + 2 * s - 1) / (4 * (s + j - 1) * (s + j))
def v(j): return -j * (j + 2 * t - 1) / (4 * (t + j - 1) * (t + j))

cS = lambda i: u(1) if i == 0 else u(2 * i) + u(2 * i + 1)
lS = lambda i: u(2 * i - 1) * u(2 * i)
cT = lambda m: v(1) if m == 0 else v(2 * m) + v(2 * m + 1)
lT = lambda m: v(2 * m - 1) * v(2 * m)

def build_OPS(c_fun, l_fun, imax):
    P = [Ry(1)]
    if imax >= 1: P.append(y - c_fun(0))
    for i in range(1, imax):
        P.append((y - c_fun(i)) * P[i] - l_fun(i) * P[i - 1])
    return P

def kappa_closed(i, m):
    d = i - m
    if d < 0 or m > i: return F(0)
    return (factorial(2 * i) / factorial(2 * m) / QQ(4) ** d
            * binomf(delta, d) / (poch(t + 2 * m + 1, d) * poch(s + i + m, d)))

def expand_in(Pi, Qlist):
    coeffs = {}; q = Pi
    while q != 0 and q.degree() >= 0:
        dg = q.degree(); lc = q.leading_coefficient()
        coeffs[dg] = lc; q = q - lc * Qlist[dg]
    return coeffs

def Sfun(p): return sum(cf * a_moment(2 * k) for k, cf in enumerate(p.list()))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    P = build_OPS(cS, lS, N + 1)
    Q = build_OPS(cT, lT, N + 1)

    print("Lemma 23.6 (lem:bessel-kappa) of d24Bessel.tex   [s,t symbolic]")
    print("P_i = sum_m kappa_{i,m} Q_m ; kappa closed form checked over QQ(s,t).")
    print("")
    print(" (i,m) |  dif = LHS(i,m) - RHS(i,m)")
    print("-------+---------------------------")

    all_zero = True
    for i in range(N + 1):
        c = expand_in(P[i], Q)
        for m in range(i + 1):
            d = F(c.get(m, F(0)) - kappa_closed(i, m))
            if d != 0: all_zero = False
            print("(%d,%d) | %s" % (i, m, d))

    sane = all(F(Sfun(P[i] * y ** j)) == 0 for i in range(1, N + 1) for j in range(i))
    print("")
    print("Orthogonality sanity S[P_i y^j]=0 (j<i):", sane)
    print("")
    if all_zero and sane:
        print("Resume: all dif = 0 for pairs (i,m), i<=%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
