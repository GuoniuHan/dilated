# check_Lem25_1.sage -- verification of Lemma 25.1 (lem:xbe-moments) of d25xBessel.tex.
#
# The multiplicative Bessel family  f_nu(x) = (1+x) cosb_nu(x)^2,  where
#   cosb_nu(x) = sum_{k>=0} (-1)^k / (4^k k! (nu+1)_k) x^{2k}.
# With a_m = m! [x^m] f_nu, the lemma asserts the even/odd moments are a single
# hypergeometric term:
#   a_{2k}   = (-4)^k (1/2)_k (nu+1/2)_k / ((nu+1)_k (2nu+1)_k),
#   a_{2k+1} = (2k+1) a_{2k} = (-4)^k (3/2)_k (nu+1/2)_k / ((nu+1)_k (2nu+1)_k).
#
# Check: build a_{2k} DIRECTLY from the definition (Cauchy square of cosb_nu, times
# (2k)!) and a_{2k+1}=(2k+1)a_{2k}, and compare to the two closed forms, as rational
# functions of the INDETERMINATE nu (dif = 0 in QQ(nu)).
#
# Run:  sage check_Lem25_1.sage N    (checks k = 0..N)

import sys

R = PolynomialRing(QQ, 'nu')
nu = R.gen()
K = R.fraction_field()

def poch(base, k):
    return prod(base + j for j in range(k))

def a2k_closed(k):
    return (K((-4) ** k) * poch(K(1) / 2, k) * poch(nu + K(1) / 2, k)
            / (poch(nu + 1, k) * poch(2 * nu + 1, k)))

def a2k1_closed(k):
    return (K((-4) ** k) * poch(K(3) / 2, k) * poch(nu + K(1) / 2, k)
            / (poch(nu + 1, k) * poch(2 * nu + 1, k)))

def a_even_from_def(K_max):
    # a_{2k} = (2k)! [x^{2k}] cosb_nu^2  by Cauchy square, over QQ(nu)
    c = [K((-1) ** k) / (K(4 ** k * factorial(k)) * poch(nu + 1, k)) for k in range(K_max + 1)]
    a = []
    for k in range(K_max + 1):
        e = sum(c[j] * c[k - j] for j in range(k + 1))
        a.append(K(factorial(2 * k)) * e)
    return a

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    adef = a_even_from_def(N)

    print("Lemma 25.1 (lem:xbe-moments) of d25xBessel.tex   [exact, nu over QQ(nu)]")
    print("")
    print("  k | dif a_{2k} (def vs closed) | dif a_{2k+1} (=(2k+1)a_{2k} vs closed)")
    print("----+----------------------------+---------------------------------------")

    all_zero = True
    for k in range(N + 1):
        d_even = adef[k] - a2k_closed(k)
        d_odd = (2 * k + 1) * adef[k] - a2k1_closed(k)
        if d_even != 0 or d_odd != 0:
            all_zero = False
        print("%3d | %26s | %s" % (k, d_even, d_odd))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for k = 0..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
