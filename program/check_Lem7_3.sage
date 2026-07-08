# check_Lem7_3.sage -- verification of Lemma 7.3 (lem:connGene) of d04euler.tex.
#
# Connection coefficients of the (s,t) family.  Let P_i be the monic
# S-orthogonal polynomials (u_j = j(j+s)) and Q_m the monic T-orthogonal
# polynomials (v_j = j(j+t)), built from the three-term recurrences with
#   c^S_i = 2i(2i+s)+(2i+1)(2i+1+s),  l^S_i = (2i-1)(2i-1+s) 2i(2i+s),
#   c^T_m = 2m(2m+t)+(2m+1)(2m+1+t),  l^T_m = (2m-1)(2m-1+t) 2m(2m+t).
# Writing P_i = sum_m kappa_{i,m} Q_m, the lemma asserts
#   kappa_{i,m} = (2i)!/(2m)! * binomial((t-s)/2, i-m),
# and consequently
#   T[P_i Q_m] = (2i)! (t+1)_{2m} binomial((t-s)/2, i-m).
#
# Two independent checks, both with s,t kept SYMBOLIC over QQ:
#   (A) kappa_{i,m} obtained by expanding P_i in the Q-basis  vs the formula;
#   (B) T[P_i Q_m] with T the moment functional mu^T_m = E_{2m}^{(t+1)}  vs formula.
#
# Indexed by pairs; "sage check_Lem7_3.sage N" checks all 0 <= i,m <= N
# (kappa_{i,m} = 0 for m > i, so the vanishing is checked too).
#
# Run:  sage check_Lem7_3.sage N
# prints, for each pair, LHS, RHS and dif = LHS-RHS,
# then a resume stating whether every dif is 0.

import sys

Rst = PolynomialRing(QQ, 's, t')
s, t = Rst.gens()
Ry = PolynomialRing(Rst, 'y')
y = Ry.gen()

def poch(x, k):
    return prod(x + j for j in range(k))

def binom_half(a, d):
    # binomial(a, d) = a(a-1)...(a-d+1)/d!  for integer d; 0 for d < 0
    if d < 0:
        return Rst(0)
    return prod(a - j for j in range(d)) / factorial(d)

# recurrence coefficients (symbolic in s, t)
cS = lambda i: 2*i*(2*i + s) + (2*i + 1)*(2*i + 1 + s)
lS = lambda i: (2*i - 1)*(2*i - 1 + s) * 2*i*(2*i + s)
cT = lambda m: 2*m*(2*m + t) + (2*m + 1)*(2*m + 1 + t)
lT = lambda m: (2*m - 1)*(2*m - 1 + t) * 2*m*(2*m + t)

def build_family(N, cfun, lfun):
    p = [Ry(1)]
    if N >= 1:
        p.append(y - cfun(0))
    for k in range(1, N):
        p.append((y - cfun(k)) * p[k] - lfun(k) * p[k - 1])
    return p

def expand_kappa(Pi, Q, i):
    # coefficients of monic P_i in the monic Q-basis:  P_i = sum_m kappa[m] Q_m
    rem = Pi
    kap = [Rst(0)] * (i + 1)
    for m in range(i, -1, -1):
        cm = rem[m]                     # y^m coeff (rem currently has degree m)
        kap[m] = cm
        rem = rem - cm * Q[m]
    return kap

def moments_T(M):
    # mu^T_m = E_{2m}^{(t+1)} = (2m)! [x^{2m}] cos(x)^{-(t+1)},  m = 0..M
    Rx = PowerSeriesRing(Rst, 'x', default_prec=2 * M + 3)
    x = Rx.gen()
    Tsec = (-(t + 1) * cos(x).log()).exp()
    return [factorial(2 * m) * Tsec[2 * m] for m in range(M + 1)]

def apply_T(poly, mu):
    return sum(coef * mu[k] for k, coef in enumerate(poly.list()))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 4

    P = build_family(N, cS, lS)
    Q = build_family(N, cT, lT)
    muT = moments_T(2 * N)
    kap_expand = [expand_kappa(P[i], Q, i) for i in range(N + 1)]

    def kappa_formula(i, m):
        return factorial(2*i) / factorial(2*m) * binom_half((t - s)/2, i - m)

    print("Lemma 7.3 (lem:connGene) of d04euler.tex   [s, t kept symbolic]")
    print("")

    all_zero = True

    print("(A) kappa_{i,m}:  expansion of P_i in Q-basis  vs  (2i)!/(2m)! C((t-s)/2, i-m)")
    print("(i,m) |                LHS                |  dif")
    print("------+-----------------------------------+-----")
    for i in range(N + 1):
        for m in range(N + 1):
            L = kap_expand[i][m] if m <= i else Rst(0)
            R = kappa_formula(i, m)
            d = Rst(L - R)
            if d != 0:
                all_zero = False
            print("(%d,%d) | %33s | %s" % (i, m, L, d))

    print("")
    print("(B) T[P_i Q_m]  vs  (2i)! (t+1)_{2m} C((t-s)/2, i-m)")
    print("(i,m) |                LHS                |  dif")
    print("------+-----------------------------------+-----")
    for i in range(N + 1):
        for m in range(N + 1):
            L = apply_T(P[i] * Q[m], muT)
            R = factorial(2*i) * poch(t + 1, 2*m) * binom_half((t - s)/2, i - m)
            d = Rst(L - R)
            if d != 0:
                all_zero = False
            print("(%d,%d) | %33s | %s" % (i, m, L, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for 0 <= i,m <= %d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
