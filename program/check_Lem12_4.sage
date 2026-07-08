# check_Lem12_4.sage -- verification of Lemma 12.4 (lem:deriv) of d06xcos.tex.
#
# Derivative connection coefficients for the monic S-orthogonal polynomials P_i
# of f(x) = (1+x)/cos(x)^{s+1}, where S[y^k] = E^{(s+1)}_{2k}.  The lemma asserts:
#
#   (A) the three-term recurrence P_{i+1} = (y - b_i) P_i - lambda_i P_{i-1} has
#         b_i     = (4i+1) s + 8i^2 + 4i + 1,
#         lambda_i = 2i(2i-1)(s+2i-1)(s+2i);
#
#   (B) with  T_i := y P_i' - i P_i = sum_{a=0}^{i-1} Delta_{i,a} P_a,
#         Delta_{i,a} = (-1)^{i-a-1} 2^{i-a-1}/(2(i-a)-1) * (i!/a!)
#                       * (2i-1)!!/(2a-1)!! * ( s + (2(i+a)+1)/(2(i-a)+1) ),
#       and in particular deg_s Delta_{i,a} <= 1.
#
# The P_i are built directly FROM the moments (Gram-Schmidt three-term recurrence),
# so (A) checks the stated b_i,lambda_i and (B) checks Delta independently.
# Everything is exact; s is an indeterminate (work over K = QQ(s)).
#
# Run:  sage check_Lem12_4.sage N
# prints check (A) for i=0..N-1 and (B) for pairs 0<=a<i<=N, then a resume.

import sys

R = PolynomialRing(QQ, 's')
s = R.gen()
K = R.fraction_field()
Ky = PolynomialRing(K, 'y')
y = Ky.gen()

def df(k):
    return prod(range(k, 0, -2))          # double factorial k!! (df(-1)=df(0)=1)

def E_moments(Kmax):
    Nser = 2 * Kmax + 2
    St = PowerSeriesRing(QQ, 'x', default_prec=Nser + 3)
    x = St.gen()
    tanser = sin(x) / cos(x)
    t = [QQ(tanser[i]) for i in range(Nser + 1)]
    g = [R(0)] * (Nser + 1)
    g[0] = R(1)
    for j in range(Nser):
        acc = R(0)
        for i in range(1, j + 1):
            acc += t[i] * g[j - i]
        g[j + 1] = (s + 1) * acc / (j + 1)
    return [factorial(2 * k) * g[2 * k] for k in range(Kmax + 1)]

def expand(P, basis):
    # coefficients c[a] of monic P in the monic basis[a] (deg basis[a]=a): P=sum c[a] basis[a]
    rem = P
    c = {}
    while rem != 0:
        dd = rem.degree()
        c[dd] = rem.leading_coefficient()
        rem = rem - c[dd] * basis[dd]
    return c

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    E = E_moments(2 * N + 2)
    muS = [K(E[t]) for t in range(len(E))]

    def Sapply(poly):
        return sum(poly[t] * muS[t] for t in range(poly.degree() + 1)) if poly != 0 else K(0)

    # monic S-orthogonal polynomials, storing recurrence coefficients
    P = [Ky(1)]
    b = []
    lam = [None]
    b0 = Sapply(y * P[0]) / Sapply(P[0])
    b.append(b0)
    P.append(y - b0)
    for i in range(1, N + 1):
        bi = Sapply(y * P[i] ** 2) / Sapply(P[i] ** 2)
        li = Sapply(P[i] ** 2) / Sapply(P[i - 1] ** 2)
        b.append(bi)
        lam.append(li)
        P.append((y - bi) * P[i] - li * P[i - 1])

    def b_formula(i):
        return K((4 * i + 1) * s + 8 * i ** 2 + 4 * i + 1)

    def lam_formula(i):
        return K(2 * i * (2 * i - 1) * (s + 2 * i - 1) * (s + 2 * i))

    def Delta_formula(i, a):
        d = i - a
        sign = (-1) ** (d - 1)
        return (K(sign) * K(2) ** (d - 1) / K(2 * d - 1)
                * K(factorial(i)) / K(factorial(a))
                * K(df(2 * i - 1)) / K(df(2 * a - 1))
                * (s + QQ(2 * (i + a) + 1) / QQ(2 * d + 1)))

    print("Lemma 12.4 (lem:deriv) of d06xcos.tex   [exact, s over QQ(s)]")
    print("")
    print("(A) recurrence:  b_i vs (4i+1)s+8i^2+4i+1 ;  lambda_i vs 2i(2i-1)(s+2i-1)(s+2i)")
    print("  i | dif b_i | dif lambda_i")
    print("----+---------+-------------")

    all_zero = True
    for i in range(0, N + 1):
        db = b[i] - b_formula(i)
        if db != 0:
            all_zero = False
        if i >= 1:
            dl = lam[i] - lam_formula(i)
            if dl != 0:
                all_zero = False
            print("%3d | %7s | %s" % (i, db, dl))
        else:
            print("%3d | %7s | %s" % (i, db, "--"))

    print("")
    print("(B) Delta_{i,a}:  expand T_i = y P_i' - i P_i in P-basis  vs closed form")
    print("(i,a) | deg_s | dif = LHS - RHS  (should be 0)")
    print("------+-------+-------------------------------")
    for i in range(1, N + 1):
        Ti = y * P[i].derivative() - i * P[i]
        coeffs = expand(Ti, P)
        for a in range(i):
            L = coeffs.get(a, K(0))
            Rn = Delta_formula(i, a)
            d = L - Rn
            if d != 0:
                all_zero = False
            degL = L.numerator().degree() if L != 0 else 0
            print("(%d,%d) | %5s | %s" % (i, a, degL, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 (A for i=0..%d, B for a<i). Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
