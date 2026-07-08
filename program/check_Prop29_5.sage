# check_Prop29_5.sage -- verification of Proposition 29.5 (prop:hankel-springer)
# of d57hankel.tex.
#
# The two-parameter Springer family.  Let a_n = n! [x^n] (cos x - t sin x)^{-r}.
# Proposition 29.5 asserts:
#   (i) the ordinary generating function sum a_n z^n has the Jacobi continued
#       fraction with  c_n = t(2n+r),  lambda_n = (1+t^2) n (n+r-1);
#   (ii) consequently, for all n>=1,
#       H_n = det( a_{i+j} )_{0<=i,j<n} = (1+t^2)^{C(n,2)} prod_{i=1}^{n-1} i! (r)_i.
# At t=r=1 the a_n are the Springer numbers S_n = 1,1,3,11,57,361,2763,... and
#   H_n = 2^{C(n,2)} prod (i!)^2 = 1,2,32,9216,84934656,...
#
# t,r are kept SYMBOLIC (indeterminates over QQ).  The moments a_n are generated
# from (cos x - t sin x)^{-r} = exp(-r log(cos x - t sin x)) as a power series
# over QQ(t,r) -- independent of the claimed J-fraction.  Part (i) is verified
# by building the monic orthogonal polynomials from the claimed c_n,lambda_n and
# checking L[y^j P_i]=0 for j<i; part (ii) checks H_n against the closed form.
#
# Run:  sage check_Prop29_5.sage N

import sys

Rtr = PolynomialRing(QQ, 't,r'); Ftr = Rtr.fraction_field()
t, r = Ftr(Rtr.gen(0)), Ftr(Rtr.gen(1))
Ry = PolynomialRing(Ftr, 'y'); y = Ry.gen()

def poch(x, k): return prod(x + j for j in range(k))

def build_a(M):
    R = PowerSeriesRing(Ftr, 'x', default_prec=M + 2)
    x = R.gen()
    base = cos(x) - t * sin(x)                 # power series, base(0)=1
    f = exp(-r * log(base))
    return [factorial(k) * f[k] for k in range(M + 1)]

def cN(n): return t * (2 * n + r)
def lamN(n): return (1 + t ** 2) * n * (n + r - 1)

def buildP(imax):
    P = [Ry(1)]
    if imax >= 1: P.append(y - cN(0))
    for n in range(1, imax):
        P.append((y - cN(n)) * P[n] - lamN(n) * P[n - 1])
    return P

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(2 * N)
    P = buildP(N + 1)

    def Lfun(p):
        return sum(cf * a[k] for k, cf in enumerate(p.list()))

    print("Proposition 29.5 (prop:hankel-springer) of d57hankel.tex   [t,r symbolic]")
    print("(i) J-fraction c_n=t(2n+r), lam_n=(1+t^2)n(n+r-1);  (ii) H_n closed form.")
    print("")

    # Part (i): orthogonality of the claimed J-fraction OPS
    orth = all(Ftr(Lfun(P[i] * y ** j)) == 0 for i in range(1, N + 1) for j in range(i))
    print("Part (i)  J-fraction orthogonality L[y^j P_i]=0 (j<i):", orth)
    print("")

    print("Part (ii)  H_n = (1+t^2)^C(n,2) prod_{i<n} i!(r)_i")
    print("  n |  dif = H_n - RHS(n)")
    print("----+--------------------")
    all_zero = orth
    for n in range(1, N + 1):
        H = matrix(Ftr, n, n, lambda i, j: a[i + j]).det()
        f1 = (1 + t ** 2) ** binomial(n, 2) * prod(factorial(i) * poch(r, i) for i in range(1, n))
        d = Ftr(H - f1)
        if d != 0:
            all_zero = False
        print("%3d | %s" % (n, d))

    # t=r=1 numeric: Springer numbers and Hankel values
    a11 = [Ftr(v).subs(t=1, r=1) for v in a]
    S_expected = [1, 1, 3, 11, 57, 361, 2763]
    s_ok = all(a11[k] == S_expected[k] for k in range(min(len(a11), len(S_expected))))
    print("")
    print("t=r=1 Springer numbers S_n match:", s_ok, " ", a11[:min(7, len(a11))])
    hv_ok = True
    for n in range(1, N + 1):
        Hn = matrix(QQ, n, n, lambda i, j: a11[i + j]).det()
        if Hn != 2 ** binomial(n, 2) * prod(factorial(i) ** 2 for i in range(1, n)):
            hv_ok = False
    print("t=r=1 Hankel values H_n = 2^C(n,2) prod(i!)^2 match:", hv_ok)
    if not (s_ok and hv_ok):
        all_zero = False

    print("")
    if all_zero:
        print("Resume: all checks pass for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some check FAILED. Formula NOT checked!")
