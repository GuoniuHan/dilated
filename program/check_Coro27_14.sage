# check_Coro27_14.sage -- Corollary 27.14 (cor:32) of d51coro.tex: Euler (s=3,t=2).
#   HH_n     = 2^{-nund} 6^{-nbar} Omega(-1) prod_i (2i)! prod_l (2l+3)! prod_m (2m+2)!
#   HH_n^(2) = (2n-1)!! 2^nbar (nbar+1)! (2nund+1)!! HH_n
# Run:  sage check_Coro27_14.sage N
import sys
def _run(title, claims, dN=6):
    N = int(sys.argv[1]) if len(sys.argv) > 1 else dN
    print(title); print(); ok = True
    for n in range(1, N+1):
        parts = []
        for nm, L, R in claims:
            d = L(n) - R(n)
            if d != 0: ok = False
            parts.append("%s (dif %s)" % (nm, d))
        print("n=%d:  %s" % (n, "   ".join(parts)))
    print(); print("Resume: all dif = 0. Formula checked." if ok else "Resume: some dif != 0. NOT checked!")
def df(k): return prod(range(k, 0, -2))
def sfrac(bfun, K):
    PS = PowerSeriesRing(QQ, 'x', default_prec=K+2); x = PS.gen(); f = PS(1)
    for j in range(K+1, 0, -1): f = (1 - bfun(j)*x*f)^(-1)
    return [f[k] for k in range(K+1)]
def Omega(delta, n):
    nb=(n+1)//2; nu=n//2
    return (-1)^binomial(nb,2)*prod((QQ(delta)/2+cc-rr)/(n-rr-cc+1) for rr in range(1,nb+1) for cc in range(1,nu+1))

s, t = 3, 2
N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 6
_MU = (sfrac(lambda j: j*(j+s), 2*N0+2), sfrac(lambda j: j*(j+t), 2*N0+2))
def a(m): return _MU[0][m//2] if m % 2 == 0 else _MU[1][m//2]
def HH(n, sh=0): return matrix(QQ, n, n, lambda i, j: a(2*i+j+sh)).det()

def base(n):
    nb=(n+1)//2; nu=n//2
    return (2^(-nu)*6^(-nb)*Omega(-1,n)*prod(factorial(2*i) for i in range(n))
            *prod(factorial(2*l+3) for l in range(nb))*prod(factorial(2*m+2) for m in range(nu)))
def sh2(n):
    nb=(n+1)//2; nu=n//2
    return df(2*n-1)*2^nb*factorial(nb+1)*df(2*nu+1)*base(n)

_run("Corollary 27.14 (cor:32): Euler (s=3,t=2)",
     [("HH_n", lambda n: HH(n), base), ("HH_n^(2)", lambda n: HH(n,2), sh2)])
