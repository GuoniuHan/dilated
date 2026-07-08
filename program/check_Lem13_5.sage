# check_Lem13_5.sage -- verification of Lemma 13.5 (lem:percol) of d07xcosshift.tex.
#
# Per-column degrees.  Ph_a are the kernel polynomials of S (the monic
# S1-orthogonal polynomials), T*[y^k] = (2k+1) E^{(s+1)}_{2k}.  The lemma asserts
#   (i)  deg_s T*[Ph_a y^l] <= 2l          for 0 <= l < a,
#   (ii) T*[Ph_a] = (-1)^a (2a)!           (the case l = 0).
#
# s is kept SYMBOLIC over QQ.  Ph_a are built as kernel polynomials
# Ph_a = (P_{a+1} - nu_a P_a)/y, nu_a = P_{a+1}(0)/P_a(0), from the monic
# S-orthogonal polynomials P_i of eq:Sdata.
#
# Run:  sage check_Lem13_5.sage N
# prints two tables:
#   (i)  for each (a,l) with 0<=l<a<=N: deg_s, the bound 2l, and OK/FAIL,
#   (ii) for each a<=N: T*[Ph_a] vs (-1)^a (2a)! with dif = LHS-RHS,
# then a resume stating whether every bound holds and every dif is 0.

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()
Fs = Rs.fraction_field()
Ry = PolynomialRing(Fs, 'y')
y = Ry.gen()

def build_moments(M):
    Rx = PowerSeriesRing(Rs, 'x', default_prec=2 * M + 4)
    x = Rx.gen()
    Ssec = (-(s + 1) * cos(x).log()).exp()
    return [factorial(2 * k) * Ssec[2 * k] for k in range(M + 1)]

cS = lambda i: (4*i + 1)*s + 8*i**2 + 4*i + 1
lS = lambda i: 2*i*(2*i - 1)*(s + 2*i - 1)*(s + 2*i)

def buildP(N):
    P = [Ry(1)]
    if N >= 1:
        P.append(y - cS(0))
    for k in range(1, N):
        P.append((y - cS(k)) * P[k] - lS(k) * P[k - 1])
    return P

def sdeg(v):
    # s-degree of an element v of Fs = QQ(s) (= numerator deg - denominator deg)
    return v.numerator().degree() - v.denominator().degree()

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    E = build_moments(2 * N + 2)
    Tsmom = lambda k: Fs((2*k + 1) * E[k])
    apply_mom = lambda p, mom: sum(c * mom(k) for k, c in enumerate(p.list()))

    P = buildP(N + 2)
    def Ph(a):
        nu = P[a + 1].subs(y=0) / P[a].subs(y=0)
        return (P[a + 1] - nu * P[a]).shift(-1)

    print("Lemma 13.5 (lem:percol) of d07xcosshift.tex   [s kept symbolic]")
    print("")

    all_ok = True

    print("(i) degree bound: deg_s T*[Ph_a y^l] <= 2l   (0 <= l < a)")
    print("(a,l) | deg_s | bound 2l | status")
    print("------+-------+----------+-------")
    for a in range(1, N + 1):
        for l in range(a):
            dg = sdeg(apply_mom(Ph(a) * y**l, Tsmom))
            ok = dg <= 2 * l
            if not ok:
                all_ok = False
            print("(%d,%d) | %5d | %8d | %s" % (a, l, dg, 2 * l, "OK" if ok else "FAIL"))

    print("")
    print("(ii) exact value: T*[Ph_a] = (-1)^a (2a)!")
    print("  a |            LHS            |            RHS            |  dif")
    print("----+--------------------------+--------------------------+-----")
    for a in range(N + 1):
        L = apply_mom(Ph(a), Tsmom)
        R = Fs((-1) ** a * factorial(2 * a))
        d = L - R
        if d != 0:
            all_ok = False
        print("%3d | %24s | %24s | %s" % (a, L, R, d))

    print("")
    if all_ok:
        print("Resume: all bounds hold and all dif = 0 for indices <= %d. Formula checked." % N)
    else:
        print("Resume: some bound fails or some dif != 0. Formula NOT checked!")
