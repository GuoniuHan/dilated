# check_Ex23_3.sage -- verification of Example 23.3 (ex:besselJ0Si) of d24Bessel.tex.
#
# At (s,t) = (0, 1/2) the Bessel (s,t) family member is
#   f_{0,1/2}(x) = cosb_0(x) + int_0^x cosb_{1/2}(y) dy = J_0(x) + Si(x),
# since cosb_0 = J_0 and int_0^x cosb_{1/2} = int_0^x sin(y)/y dy = Si(x).
# Here 2(t-s)=1, so the signed factor is Omega(1)=2^{-C(n,2)}, and
# Theorem 23.1 specialises to the factorial product
#   HH_n(J0+Si) = (prod_{i<n}(2i)!)/(2^{n(n-1)} 2^{C(n,2)})
#                 * prod_{l=0}^{nb-1} (1/2)_l/(n-1+l)!
#                 * prod_{m=0}^{nn-1} m!/(3/2)_{n-1+m}
#               = 1, 1/6, 1/960, 1/1693440, ...   (n=1,2,3,4).
#
# The moments are exact rationals at (s,t)=(0,1/2).  Each HH_n is checked in
# three ways: (i) the direct determinant; (ii) the tabulated values; (iii) the
# specialisation of Theorem 23.1's closed form at (0,1/2).
#
# Run:  sage check_Ex23_3.sage N

import sys

H = QQ(1) / 2

def poch(x, k): return prod(x + j for j in range(k))

def a_moment(midx):
    # (s,t) = (0, 1/2)
    if midx % 2 == 0:
        k = midx // 2
        return (-1) ** k * poch(H, k) / poch(QQ(0) + 1, k)
    k = (midx - 1) // 2
    return (-1) ** k * poch(H, k) / poch(H + 1, k)

def HH(n):
    return matrix(QQ, n, n, lambda i, j: a_moment(2 * i + j)).det()

def prodform(n):
    nb = (n + 1) // 2; nn = n // 2
    val = prod(factorial(2 * i) for i in range(n)) / (QQ(2) ** (n * (n - 1)) * QQ(2) ** binomial(n, 2))
    val *= prod(poch(H, l) / factorial(n - 1 + l) for l in range(nb))
    val *= prod(factorial(m) / poch(QQ(3) / 2, n - 1 + m) for m in range(nn))
    return val

def thm_closed_at_0_half(n):
    # Theorem 23.1 RHS specialised at (s,t)=(0,1/2)
    nb = (n + 1) // 2; nn = n // 2
    s0 = QQ(0); t0 = H
    val = (-1) ** binomial(nb, 2) * prod(factorial(2 * i) for i in range(n)) / QQ(2) ** (n * (n - 1))
    val *= prod(poch(s0 + H, l) / poch(s0 + 1, n - 1 + l) for l in range(nb))
    val *= prod(poch(t0 + H, m) / poch(t0 + 1, n - 1 + m) for m in range(nn))
    val *= prod((t0 - s0 + c - r) / (n - r - c + 1)
                for r in range(1, nb + 1) for c in range(1, nn + 1))
    return val

VALS = [1, QQ(1) / 6, QQ(1) / 960, QQ(1) / 1693440]

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 4

    print("Example 23.3 (ex:besselJ0Si) of d24Bessel.tex   [(s,t)=(0,1/2), J0+Si]")
    print("")
    print("  n |     HH_n(J0+Si)     | dif(list) | dif(prodform) | dif(thm@0,1/2)")
    print("----+---------------------+-----------+---------------+---------------")

    all_zero = True
    for n in range(1, N + 1):
        L = HH(n)
        d_list = (L - VALS[n - 1]) if n - 1 < len(VALS) else QQ(0)
        d_pf = L - prodform(n)
        d_thm = L - thm_closed_at_0_half(n)
        if d_list != 0 or d_pf != 0 or d_thm != 0:
            all_zero = False
        print("%3d | %19s | %9s | %13s | %s" % (n, L, d_list, d_pf, d_thm))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
