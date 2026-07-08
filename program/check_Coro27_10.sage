# check_Coro27_10.sage -- Corollary 27.10 (cor:12) of d51coro.tex: Euler (s=1,t=2).
#   HH_n     = 2^{-C(n,2)-nund} prod_{k<n}(2k)! prod_{j=1}^{n} j!
#   HH_n^(2) = (2n-1)!! 2^nbar nbar! (2nund+1)!! HH_n
# Run:  sage check_Coro27_10.sage N
import sys
def _run(title, claims, dN=6):
    N = int(sys.argv[1]) if len(sys.argv) > 1 else dN
    print(title); print(); ok = True
    for n in range(1, N+1):
        parts = [ ("%s (dif %s)" % (nm, (L(n)-R(n)))) for nm, L, R in claims]
        if any((L(n)-R(n)) != 0 for nm,L,R in claims): ok = False
        print("n=%d:  %s" % (n, "   ".join(parts)))
    print(); print("Resume: all dif = 0. Formula checked." if ok else "Resume: some dif != 0. NOT checked!")

s, t = 1, 2
def df(k): return prod(range(k, 0, -2))
def sfrac(bfun, K):
    PS = PowerSeriesRing(QQ, 'x', default_prec=K+2); x = PS.gen(); f = PS(1)
    for j in range(K+1, 0, -1): f = (1 - bfun(j)*x*f)^(-1)
    return [f[k] for k in range(K+1)]
N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 6
_MU = (sfrac(lambda j: j*(j+s), 2*N0+2), sfrac(lambda j: j*(j+t), 2*N0+2))
def a(m): return _MU[0][m//2] if m % 2 == 0 else _MU[1][m//2]
def HH(n, sh=0): return matrix(QQ, n, n, lambda i, j: a(2*i+j+sh)).det()

def base(n):
    nu=n//2
    return 2^(-binomial(n,2)-nu)*prod(factorial(2*k) for k in range(1,n))*prod(factorial(j) for j in range(1,n+1))
def sh2(n):
    nb=(n+1)//2; nu=n//2
    return df(2*n-1)*2^nb*factorial(nb)*df(2*nu+1)*base(n)

_run("Corollary 27.10 (cor:12): Euler (s=1,t=2)",
     [("HH_n", lambda n: HH(n), base), ("HH_n^(2)", lambda n: HH(n,2), sh2)])
