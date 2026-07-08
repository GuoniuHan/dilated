# check_Lem12_2.sage -- verification of Lemma 12.2 (lem:reduction) of d06xcos.tex.
#
# One-sided reduction for f(x) = (1+x)/cos(x)^{s+1}, general s.  With
#   S[y^k]  = E^{(s+1)}_{2k}        (even functional, classical),
#   T*[y^k] = (2k+1) E^{(s+1)}_{2k} (odd functional),
# let P_0,P_1,... be the monic S-orthogonal polynomials and
#   h^S_l = S[P_l^2] = (2l)! (s+1)_{2l}.
# The lemma asserts (bar n = ceil(n/2), un n = floor(n/2)):
#
#   (A)  HH_n(f) = (-1)^C(bar n,2) ( prod_{l=0}^{bar n-1} h^S_l )
#                  * det( T*[ P_{bar n + r} y^m ] )_{0<=r,m<=un n-1}.        [eq:reduction]
#
#   (B)  T* = S o (2 y d/dy + 1), so for i > m
#           T*[ P_i y^m ] = 2 S[ y^{m+1} P_i' ]                              [eq:Dentry]
#        (the term (2m+1) S[y^m P_i] vanishes by orthogonality since deg y^m < i).
#
# Everything is exact; s is kept as an indeterminate (work over K = QQ(s)).
#
# Run:  sage check_Lem12_2.sage N
# prints checks (A) for n=1..N and (B) for pairs i>m, then a resume.

import sys

R = PolynomialRing(QQ, 's')
s = R.gen()
K = R.fraction_field()
Ky = PolynomialRing(K, 'y')
y = Ky.gen()

def df(k):
    return prod(range(k, 0, -2))

def E_moments(Kmax):
    # E^{(s+1)}_{2k} = (2k)! [x^{2k}] sec(x)^{s+1}, k=0..Kmax, over QQ[s].
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

def poch(base, i):
    return prod(base + j for j in range(i))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    E = E_moments(2 * N + 2)
    muS = [K(E[t]) for t in range(len(E))]                 # S[y^t]
    muT = [K((2 * t + 1) * E[t]) for t in range(len(E))]   # T*[y^t]

    def Sapply(poly):
        return sum(poly[t] * muS[t] for t in range(poly.degree() + 1)) if poly != 0 else K(0)

    def Tapply(poly):
        return sum(poly[t] * muT[t] for t in range(poly.degree() + 1)) if poly != 0 else K(0)

    # monic S-orthogonal polynomials P_0..P_N
    P = [Ky(1)]
    b0 = Sapply(y * P[0]) / Sapply(P[0])
    P.append(y - b0)
    for i in range(1, N):
        bi = Sapply(y * P[i] ** 2) / Sapply(P[i] ** 2)
        li = Sapply(P[i] ** 2) / Sapply(P[i - 1] ** 2)
        P.append((y - bi) * P[i] - li * P[i - 1])

    def hS(l):
        return K(factorial(2 * l) * poch(s + 1, 2 * l))

    def HH(n):
        def a(m):
            k = m // 2
            return muS[k] if m % 2 == 0 else muT[k]
        return matrix(K, n, n, lambda i, j: a(2 * i + j)).det()

    print("Lemma 12.2 (lem:reduction) of d06xcos.tex   [exact, s over QQ(s)]")
    print("")
    print("(A) reduction:  HH_n  vs  (-1)^C(bar n,2) prod h^S_l * det( T*[P_{bar n+r} y^m] )")
    print("  n | dif = LHS - RHS  (should be 0)")
    print("----+--------------------------------")

    all_zero = True
    for n in range(1, N + 1):
        bn = (n + 1) // 2      # ceil(n/2)
        un = n // 2            # floor(n/2)
        L = HH(n)
        if un == 0:
            detD = K(1)
        else:
            D = matrix(K, un, un,
                       lambda r, m: Tapply(P[bn + r] * y ** m))
            detD = D.det()
        Rn = K(-1) ** binomial(bn, 2) * prod(hS(l) for l in range(bn)) * detD
        d = L - Rn
        if d != 0:
            all_zero = False
        print("%3d | %s" % (n, d))

    print("")
    print("(B) T*[P_i y^m] = 2 S[y^{m+1} P_i']   for i > m")
    print("(i,m) |          LHS          |          RHS          | dif")
    print("------+-----------------------+-----------------------+----")
    for i in range(1, N + 1):
        for m in range(i):          # m < i
            L = Tapply(P[i] * y ** m)
            Rn = 2 * Sapply(y ** (m + 1) * P[i].derivative())
            d = L - Rn
            if d != 0:
                all_zero = False
            print("(%d,%d) | %21s | %21s | %s" % (i, m, L, Rn, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 (A for n=1..%d, B for i>m). Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
