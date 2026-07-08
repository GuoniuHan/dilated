# check_Thm18_1.sage -- verification of Theorem 18.1 (conj:xsin-closed) of d12xsinx.tex.
#
# Reciprocal-sine case f(x) = (1+x) x/sin x.  With g = x/sin x,
#   b_k = (2k)! [x^{2k}] g,   a_{2k} = b_k,   a_{2k+1} = (2k+1) b_k,
# and HH_n = det(a_{2i+j}).  Theorem 18.1: with
#   Q(k) = 144*64^k (k!)^5 ((k+1)!)^6 (2k)!(2k+1)! / ((3k)!(3k+3)!(3k+4)!),
# the determinants satisfy HH_{n+2} = (2/3) Q(n) HH_n, hence
#   HH_{2N}   = (2/3)^N prod_{k=0}^{N-1} Q(2k),
#   HH_{2N+1} = (2/3)^N prod_{k=0}^{N-1} Q(2k+1).
#
# LHS(n) = det(a_{2i+j}); RHS(n) = the closed product form; also the two-step
# recurrence is checked (column dif2).  Exact rational arithmetic.
#
# Run:  sage check_Thm18_1.sage N

import sys

def bmoments(Kmax):
    prec = 2*Kmax + 2
    PS = PowerSeriesRing(QQ, 'x', default_prec=prec + 1)
    x = PS.gen()
    sinc = sum((-1)^k * x^(2*k) / factorial(2*k+1) for k in range(prec//2 + 1))  # sin x / x
    g = sinc.inverse()                                                          # x / sin x
    return [factorial(2*k) * g[2*k] for k in range(Kmax + 1)]

def Q(k):
    return (144 * 64^k * factorial(k)^5 * factorial(k+1)^6 * factorial(2*k) * factorial(2*k+1)
            / (factorial(3*k) * factorial(3*k+3) * factorial(3*k+4)))

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 6

    b = bmoments(3*(N+2))
    def a(m): return b[m//2] if m % 2 == 0 else (m)*b[m//2]

    def HH(n):
        return matrix(QQ, n, n, lambda i, j: a(2*i + j)).det() if n >= 1 else QQ(1)

    def RHSclosed(n):
        M = n//2
        return (QQ(2)/3)^M * prod(Q(2*k + (n % 2)) for k in range(M))

    print("Theorem 18.1 (conj:xsin-closed) of d12xsinx.tex")
    print("")
    print("  n |              LHS(n)=HH_n              | dif=LHS-closed | dif2=HH_{n+2}-(2/3)Q(n)HH_n")
    print("----+--------------------------------------+----------------+----------------------------")

    all_zero = True
    for n in range(1, N + 1):
        L = HH(n)
        d = L - RHSclosed(n)
        d2 = HH(n+2) - (QQ(2)/3)*Q(n)*HH(n)
        if d != 0 or d2 != 0: all_zero = False
        print("%3d | %36s | %14s | %s" % (n, L, d, d2))

    print("")
    if all_zero:
        print("Resume: all dif = 0 for n = 1..%d. Formula checked." % N)
    else:
        print("Resume: some dif != 0. Formula NOT checked!")
