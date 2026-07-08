# check_Coro27_20.sage -- Corollary 27.20 (cor:xc3) of d51coro.tex: secant s=3, f=(1+x)/cos^4 x.
#   HH_n     = c_n prod_{i=1}^{n-1}(i+3)!/6
#   HH_n^(1) = ((n-1)!!)^2 HH_n
# Run:  sage check_Coro27_20.sage N
import sys
def _run(title, claims, dN=6):
    N = int(sys.argv[1]) if len(sys.argv) > 1 else dN
    print(title); print(); ok = True
    for n in range(1, N+1):
        parts = []
        for nm, L, R in claims:
            Lv = L(n); d = Lv - R(n)
            if d != 0: ok = False
            parts.append("%s=%s (dif %s)" % (nm, Lv, d))
        print("n=%d:  %s" % (n, "   ".join(parts)))
    print(); print("Resume: all dif = 0. Formula checked." if ok else "Resume: some dif != 0. NOT checked!")
def df(k): return prod(range(k, 0, -2))
def secmoments(spow, Mmax):
    prec = Mmax + 1
    PS = PowerSeriesRing(QQ, 'x', default_prec=prec+1); x = PS.gen()
    cosx = sum((-1)^k*x^(2*k)/factorial(2*k) for k in range(prec//2+1))
    f = (1+x)*cosx.inverse()^spow
    return [factorial(m)*f[m] for m in range(Mmax+1)]

N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 6
_A = secmoments(4, 3*N0+2)                      # s+1 = 4
def a(m): return _A[m]
def HH(n, sh=0): return matrix(QQ, n, n, lambda i, j: a(2*i+j+sh)).det()
def cn(n): return (2^binomial(n,2)*df(n-1)^2*prod(df(k)^6 for k in range(1, n-1))
                   / prod(factorial(i) for i in range(1, n)))
def base(n): return cn(n)*prod(factorial(i+3)/6 for i in range(1, n))
def sh1(n):  return df(n-1)^2*base(n)

_run("Corollary 27.20 (cor:xc3): secant s=3, f=(1+x)/cos^4 x",
     [("HH_n", lambda n: HH(n), base), ("HH_n^(1)", lambda n: HH(n,1), sh1)])
