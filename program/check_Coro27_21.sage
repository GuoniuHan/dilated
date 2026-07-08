# check_Coro27_21.sage -- Corollary 27.21 (rank-one perturbation, s=-2) of d51coro.tex.
# f_s = (sin x + 1)/cos^2 x + s sin x;  HH_n(f_s) = (1 - s C(n,2)) prod_{k<n}(k!)^2 (2k+1)!!.
# Run:  sage check_Coro27_21.sage N
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
sval = -2
def moments(Mmax):
    prec = Mmax + 1
    PS = PowerSeriesRing(QQ, 'x', default_prec=prec+1); x = PS.gen()
    cosx = sum((-1)^k*x^(2*k)/factorial(2*k) for k in range(prec//2+1))
    sinx = sum((-1)^k*x^(2*k+1)/factorial(2*k+1) for k in range(prec//2+1))
    f = (sinx+1)*cosx.inverse()^2 + sval*sinx
    return [factorial(m)*f[m] for m in range(Mmax+1)]
N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 6
_A = moments(3*N0+2)
def a(m): return _A[m]
def HH(n): return matrix(QQ, n, n, lambda i, j: a(2*i+j)).det()
def HE(n): return prod(factorial(k)^2*df(2*k+1) for k in range(1, n))
def base(n): return (n^2-n+1)*HE(n)
_run("Corollary 27.21 (rank-one, s=-2): f_s=(sin x+1)/cos^2 x + s sin x",
     [("HH_n", lambda n: HH(n), base)])
