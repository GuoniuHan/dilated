# check_Lem7_2.sage -- verification of Lemma 7.2 (lem:cf) of d04euler.tex.
#
# Secant and tangent polynomials (the family at s = 0, t = 1).
#
# Let S be the moment functional in the variable y with
#   S[y^k] = mu^S_k = E_{2k}   (even Euler / secant numbers 1,1,5,61,...),
# and T the moment functional with
#   T[y^m] = mu^T_m = E_{2m+1} (odd Euler / tangent numbers 1,2,16,272,...).
# The monic S-orthogonal polynomials Ph_i (= P_i|_{s=0}) and the monic
# T-orthogonal polynomials r_m (= Q_m|_{t=1}) are built from the recurrences
#   Ph_{i+1} = (y - (2i)^2 - (2i+1)^2) Ph_i - ((2i-1)(2i))^2 Ph_{i-1},
#   r_{m+1}  = (y - 2(2m+1)^2) r_m - (2m-1)(2m)^2(2m+1) r_{m-1}.
# The lemma asserts the orthogonality relations
#   S[Ph_i Ph_l] = delta_{il} ((2i)!)^2,   T[r_m r_{m'}] = delta_{mm'} (2m)!(2m+1)!.
#
# Everything is exact over QQ.  Indexed by pairs; "sage check_Lem7_2.sage N"
# checks all 0 <= i,l <= N (family S) and all 0 <= m,m' <= N (family T).
#
# Run:  sage check_Lem7_2.sage N
# prints, for each pair, LHS, RHS and dif = LHS-RHS,
# then a resume stating whether every dif is 0.

import sys

def moments(M):
    # mu^S_k = E_{2k} (k=0..M) and mu^T_m = E_{2m+1} (m=0..M), exact over QQ.
    Rx = PowerSeriesRing(QQ, 'x', default_prec=2 * M + 3)
    x = Rx.gen()
    sec  = 1 / cos(x)          # sum E_{2k} x^{2k}/(2k)!
    sec2 = sec ** 2            # (tan)' = sum E_{2k+1} x^{2k}/(2k)!
    muS = [factorial(2 * k) * sec[2 * k]  for k in range(M + 1)]
    muT = [factorial(2 * m) * sec2[2 * m] for m in range(M + 1)]
    return muS, muT

Ry = PolynomialRing(QQ, 'y')
y = Ry.gen()

def build_family(N, cfun, lfun):
    # monic orthogonal polynomials p_0..p_N via p_{k+1} = (y-c_k)p_k - l_k p_{k-1}
    p = [Ry(1)]
    if N >= 1:
        p.append(y - cfun(0))
    for k in range(1, N):
        p.append((y - cfun(k)) * p[k] - lfun(k) * p[k - 1])
    return p

def apply_functional(poly, mu):
    # apply the moment functional with moments mu to a polynomial in y
    return sum(coef * mu[k] for k, coef in enumerate(poly.list()))

# S-family (secant, s=0):  c^S_i = (2i)^2+(2i+1)^2,  l^S_i = ((2i-1)(2i))^2
cS = lambda i: (2*i)**2 + (2*i + 1)**2
lS = lambda i: ((2*i - 1) * (2*i))**2
# T-family (tangent, t=1): c^T_m = 2(2m+1)^2,  l^T_m = (2m-1)(2m)^2(2m+1)
cT = lambda m: 2 * (2*m + 1)**2
lT = lambda m: (2*m - 1) * (2*m)**2 * (2*m + 1)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    muS, muT = moments(2 * N)
    Ph = build_family(N, cS, lS)
    rr = build_family(N, cT, lT)

    print("Lemma 7.2 (lem:cf) of d04euler.tex   [exact over QQ]")
    print("secant numbers  mu^S_k = %s ..." % muS[:min(N + 1, 5)])
    print("tangent numbers mu^T_m = %s ..." % muT[:min(N + 1, 5)])
    print("")

    all_zero = True

    print("Family S:  S[Ph_i Ph_l]  vs  delta_{il} ((2i)!)^2")
    print("(i,l) |          LHS          |          RHS          |  dif")
    print("------+-----------------------+-----------------------+-----")
    for i in range(N + 1):
        for l in range(N + 1):
            L = apply_functional(Ph[i] * Ph[l], muS)
            R = factorial(2 * i)**2 if i == l else 0
            d = L - R
            if d != 0:
                all_zero = False
            print("(%d,%d) | %21s | %21s | %s" % (i, l, L, R, d))

    print("")
    print("Family T:  T[r_m r_m']  vs  delta_{mm'} (2m)!(2m+1)!")
    print("(m,m')| %21s | %21s | dif" % ("LHS", "RHS"))
    print("------+-----------------------+-----------------------+-----")
    for m in range(N + 1):
        for mp in range(N + 1):
            L = apply_functional(rr[m] * rr[mp], muT)
            R = factorial(2 * m) * factorial(2 * m + 1) if m == mp else 0
            d = L - R
            if d != 0:
                all_zero = False
            print("(%d,%d) | %21s | %21s | %s" % (m, mp, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= indices <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
