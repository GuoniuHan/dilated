# check_Lem15_4.sage -- verification of Lemma 15.4 (lem:evalminus1) of d10runkone.tex.
#
# Evaluations at v = -1.  The tangent polynomials r_m (Lemma lem:cf, written in v)
# obey  r_{m+1} = (v - cT_m) r_m - lamT_m r_{m-1},  cT_m = 2(2m+1)^2,
#   lamT_m = (2m-1)(2m)^2(2m+1);  and the shifted-secant polynomials P^(1)_l
# (Lemma 15.3) obey  P^(1)_{l+1} = (v - (2l+1)^2 - (2l+2)^2) P^(1)_l - ((2l)(2l+1))^2 P^(1)_{l-1}.
# The lemma asserts
#   r_m(-1)      = (-1)^m (2m-1)!! (2m+1)!!,
#   P^(1)_l(-1)  = (-1)^l (2l+1)!,
#   and in particular  r_m(-1)/(2m+1)! = binomial(-1/2, m).
#
# We build both families from their recurrences, evaluate at v=-1, and compare to
# the closed forms.  Everything is exact over QQ.  Indexed by m (resp. l) = 0..N.
#
# Run:  sage check_Lem15_4.sage N

import sys

def df(k):
    return prod(range(k, 0, -2))          # double factorial k!!

def gen_binom(a, k):
    # binomial(a, k) = a(a-1)...(a-k+1)/k! for k >= 0; 0 for k < 0
    if k < 0:
        return QQ(0)
    return prod(a - j for j in range(k)) / factorial(k)

Rv = PolynomialRing(QQ, 'v')
v = Rv.gen()

# tangent polynomials r_m
cT = lambda m: 2 * (2 * m + 1) ** 2
lT = lambda m: (2 * m - 1) * (2 * m) ** 2 * (2 * m + 1)
# shifted-secant polynomials P^(1)_l
cP = lambda l: (2 * l + 1) ** 2 + (2 * l + 2) ** 2
lP = lambda l: ((2 * l) * (2 * l + 1)) ** 2

def build_family(N, cfun, lfun):
    p = [Rv(1)]
    if N >= 1:
        p.append(v - cfun(0))
    for k in range(1, N):
        p.append((v - cfun(k)) * p[k] - lfun(k) * p[k - 1])
    return p

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    r = build_family(N, cT, lT)
    P = build_family(N, cP, lP)

    print("Lemma 15.4 (lem:evalminus1) of d10runkone.tex   [exact over QQ]")
    print("")

    all_zero = True

    print("(A) r_m(-1)  vs  (-1)^m (2m-1)!!(2m+1)!!   [and r_m(-1)/(2m+1)! = C(-1/2,m)]")
    print("  m |          LHS          |          RHS          | dif | binom-check")
    print("----+-----------------------+-----------------------+-----+------------")
    for m in range(N + 1):
        L = r[m](-1)
        Rn = QQ((-1) ** m) * df(2 * m - 1) * df(2 * m + 1)
        d = L - Rn
        bcheck = (L / factorial(2 * m + 1) == gen_binom(QQ(-1) / 2, m))
        if d != 0 or not bcheck:
            all_zero = False
        print("%3d | %21s | %21s | %s | %s" % (m, L, Rn, d, bcheck))

    print("")
    print("(B) P^(1)_l(-1)  vs  (-1)^l (2l+1)!")
    print("  l |          LHS          |          RHS          | dif")
    print("----+-----------------------+-----------------------+----")
    for l in range(N + 1):
        L = P[l](-1)
        Rn = QQ((-1) ** l) * factorial(2 * l + 1)
        d = L - Rn
        if d != 0:
            all_zero = False
        print("%3d | %21s | %21s | %s" % (l, L, Rn, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= m,l <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
