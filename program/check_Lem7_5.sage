# check_Lem7_5.sage -- verification of Lemma 7.5 (lem:Pahalf) of d04euler.tex.
#
# The signed factor Omega(delta), defined for nb = ceil(n/2), nn = floor(n/2) by
#   Omega(delta) = (-1)^C(nb,2) prod_{r=1}^{nb} prod_{c=1}^{nn}
#                   (delta/2 + c - r)/(n - r - c + 1).
# The double factorial is extended to negative odd arguments by
#   (2k-1)!! = (2k+1)!!/(2k+1),   so (-1)!! = 1, (-3)!! = -1, (-5)!! = 1/3, ...
#
# Two checks, both exact over QQ:
#   (A) [eq:Pa:dfact]  for every ODD delta,
#        Omega(delta) = (-1)^C(nb,2) 2^{-nb*nn} / prod_{j=nb}^{n-1} j!
#                        * prod_{j=0}^{nn-1} j! (delta+2j)!! / (delta-2nb+2j)!!.
#   (B) [eq:Pa:even / eq:Pa:odd]  the closed values for delta in {1,-1,3,-3}.
#
# In each case LHS = Omega(delta) from the defining double product; RHS = the
# claimed closed form.  Indexed by pairs (n, delta).
#
# Run:  sage check_Lem7_5.sage N
# prints, for each pair, LHS, RHS and dif = LHS-RHS,
# then a resume stating whether every dif is 0.

import sys

def dfacto(w):
    # double factorial of an ODD integer w, with the negative-odd convention
    #   w!! = (w+2)!! / (w+2)   for w <= -3.
    if w >= -1:
        r = QQ(1)
        x = w
        while x > 0:
            r *= x
            x -= 2
        return r
    return dfacto(w + 2) / QQ(w + 2)

def nb_nn(n):
    return (n + 1) // 2, n // 2          # ceil(n/2), floor(n/2)

def Omega_def(n, delta):
    # defining double product
    nb, nn = nb_nn(n)
    sign = (-1) ** binomial(nb, 2)
    prod_rc = prod((QQ(delta) / 2 + c - r) / QQ(n - r - c + 1)
                   for r in range(1, nb + 1) for c in range(1, nn + 1))
    return sign * prod_rc

def Omega_dfact(n, delta):
    # eq:Pa:dfact  (delta odd)
    nb, nn = nb_nn(n)
    sign = (-1) ** binomial(nb, 2)
    pw   = QQ(2) ** (-nb * nn)
    denom_fac = prod(factorial(j) for j in range(nb, n))
    body = prod(factorial(j) * dfacto(delta + 2*j) / dfacto(delta - 2*nb + 2*j)
                for j in range(nn))
    return sign * pw / denom_fac * body

def Omega_special(n, delta):
    # eq:Pa:even (n even) / eq:Pa:odd (n odd), for delta in {1,-1,3,-3}
    base = QQ(2) ** (-binomial(n, 2))
    if n % 2 == 0:
        m = n // 2
        return {1:  base,
                -1: (-1)**m * base,
                3:  (-1)**(m-1) * (n**2 - 1) * base,
                -3: -(n**2 - 1) * base}[delta]
    else:
        m = (n - 1) // 2
        return {1:  base,
                -1: (-1)**m * n * base,
                3:  (-1)**m * n * base,
                -3: -QQ(1)/3 * n**2 * (n**2 - 4) * base}[delta]

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    ODD = [-5, -3, -1, 1, 3, 5]          # sample odd deltas for check (A)

    print("Lemma 7.5 (lem:Pahalf) of d04euler.tex   [exact over QQ]")
    print("")

    all_zero = True

    print("(A) eq:Pa:dfact -- double-factorial formula, delta odd")
    print("(n,delta) |          Omega(def)         |          Omega(dfact)       |  dif")
    print("----------+-----------------------------+-----------------------------+-----")
    for n in range(1, N + 1):
        for delta in ODD:
            L = Omega_def(n, delta)
            R = Omega_dfact(n, delta)
            d = L - R
            if d != 0:
                all_zero = False
            print("(%2d,%3d) | %27s | %27s | %s" % (n, delta, L, R, d))

    print("")
    print("(B) eq:Pa:even / eq:Pa:odd -- special values, delta in {1,-1,3,-3}")
    print("(n,delta) |          Omega(def)         |          Omega(special)     |  dif")
    print("----------+-----------------------------+-----------------------------+-----")
    for n in range(1, N + 1):
        for delta in (1, -1, 3, -3):
            L = Omega_def(n, delta)
            R = Omega_special(n, delta)
            d = L - R
            if d != 0:
                all_zero = False
            print("(%2d,%3d) | %27s | %27s | %s" % (n, delta, L, R, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
