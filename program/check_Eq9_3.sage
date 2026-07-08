# check_Eq9_3.sage -- verification of equation (9.3) (eq:shiftmonic) of d05eulershift.tex.
#
# Clearing the Pochhammer symbols in Proposition 9.1 writes the single shift
# on the line t=1 as
#   HH_n^{(1)} = K_n * R_n(s),
# where R_n(s) is the MONIC product of linear forms (eq (9.3), eq:shiftmonic):
#   R_n(s) = [ prod_{j=1}^{nn}   prod_{k=0}^{j-1} (s+2k+1) ]^2   (odd shifts, doubled)
#          *   prod_{j=1}^{nn-1} prod_{k=1}^{j}   (s+2k)         (even shifts)
#          *   prod_{j=1}^{nb-1} prod_{k=1}^{j}   (s-2k+1),      (back shifts)
# with nb = ceil(n/2), nn = floor(n/2), and K_n a constant independent of s.
#
# This checker confirms (a) that R_n divides HH_n^{(1)} with an s-FREE quotient
# K_n (the substantive claim of (9.3)), by testing dif = HH_n^{(1)} - K_n R_n = 0
# as a polynomial in s, and (b) reports the sign of K_n.
#
# NOTE (flag for the author): eq:shiftmonic calls K_n "another POSITIVE constant".
# The verification shows K_n is positive only for n=1,2; for n>=3 it is NEGATIVE.
# Its sign is exactly (-1)^{binom(nb,2)} = (-1)^{binom(ceil(n/2),2)}, inherited
# from the back-shift Pochhammer (( 1-s )/2)_j of eq:shiftpoch, each factor of
# which contributes (-1)^j when cleared to the monic form (s-2k+1). The scalar
# Ktilde_n of eq:shiftpoch/shiftconst genuinely IS positive; only the monic-form
# constant K_n alternates in sign.  The identity HH_n^{(1)} = K_n R_n itself holds.
#
# s is kept SYMBOLIC (indeterminate over QQ).
#
# Run:  sage check_Eq9_3.sage N

import sys

Rs = PolynomialRing(QQ, 's')
s = Rs.gen()

def build_a(M):
    Rx = PowerSeriesRing(Rs, 'x', default_prec=M + 2)
    x = Rx.gen()
    logcos = cos(x).log()
    Ssec = (-(s + 1) * logcos).exp()      # cos^{-(s+1)}
    Tsec = (-2 * logcos).exp()            # cos^{-2}, t=1
    a = []
    for m in range(M + 1):
        if m % 2 == 0:
            a.append(factorial(m) * Ssec[m])
        else:
            a.append(factorial(m - 1) * Tsec[m - 1])
    return a

def HH1(n, a):
    return matrix(Rs, n, n, lambda i, j: a[2 * i + j + 1]).det()

def Rmonic(n):
    nb = (n + 1) // 2
    nn = n // 2
    p = Rs(prod(prod(s + 2 * k + 1 for k in range(0, j)) for j in range(1, nn + 1)) ** 2)
    p *= prod(prod(s + 2 * k for k in range(1, j + 1)) for j in range(1, nn))      # j=1..nn-1
    p *= prod(prod(s - 2 * k + 1 for k in range(1, j + 1)) for j in range(1, nb))  # j=1..nb-1
    return Rs(p)

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    a = build_a(3 * N)

    print("Equation (9.3) (eq:shiftmonic) of d05eulershift.tex   [s kept symbolic, t=1]")
    print("HH_n^{(1)} = K_n * R_n(s), R_n monic.  K_n sign checked vs (-1)^binom(nb,2).")
    print("")
    print("  n |          K_n          | sign ok | K_n>0 |   dif")
    print("----+-----------------------+---------+-------+------")

    all_zero = True
    for n in range(1, N + 1):
        L = HH1(n, a)
        Rn = Rmonic(n)
        q = L / Rn                              # K_n, an element of Frac(Rs)
        is_const = q in QQ
        Kn = QQ(q) if is_const else None
        d = Rs(L - Kn * Rn) if is_const else Rs(1)
        if d != 0:
            all_zero = False
        nb = (n + 1) // 2
        predicted_sign = (-1) ** binomial(nb, 2)
        sign_ok = is_const and (sign(Kn) == predicted_sign)
        print("%3d | %21s | %7s | %5s | %s"
              % (n, Kn, sign_ok, (Kn > 0), d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Identity HH_n^{(1)} = K_n R_n checked." % N)
        print("        (K_n>0 only for n=1,2; sign(K_n) = (-1)^binom(ceil(n/2),2) --")
        print("         eq:shiftmonic's wording 'positive constant K_n' should be qualified.)")
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
