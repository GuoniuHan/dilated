# check_Coro27_15.sage -- Corollary 27.15 (cor:03) of d51coro.tex: Euler (s=0,t=3).
#   HH_n     = 6^{-nund} Omega(3) prod_i (2i)! prod_l (2l)! prod_m (2m+3)!
#   HH_n^(1) = Lambda_n(0) Gamma_n(0),  Gamma_n(0)=n(n+1)-1-4nbar   (Lambda_n(0)=eq:sHH03)
#   HH_n^(2) = (2n-1)!!(2nbar-1)!! 2^nund (nund+1)! HH_n
# Run:  sage check_Coro27_15.sage N
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

s, t = 0, 3
N0 = int(sys.argv[1]) if len(sys.argv) > 1 else 6
_MU = (sfrac(lambda j: j*(j+s), 2*N0+2), sfrac(lambda j: j*(j+t), 2*N0+2))
def a(m): return _MU[0][m//2] if m % 2 == 0 else _MU[1][m//2]
def HH(n, sh=0): return matrix(QQ, n, n, lambda i, j: a(2*i+j+sh)).det()

def base(n):
    nb=(n+1)//2; nu=n//2
    return (6^(-nu)*Omega(3,n)*prod(factorial(2*i) for i in range(n))
            *prod(factorial(2*l) for l in range(nb))*prod(factorial(2*m+3) for m in range(nu)))
def Lam0(n):
    nb=(n+1)//2; nu=n//2
    rho = QQ(1)/(n+1) if n % 2 == 0 else QQ(n)/(n+2)
    return ((-1)^(nu+1)*rho*2^(-binomial(n,2))*6^(-nb)
            *prod(factorial(2*l)*factorial(2*l+3) for l in range(nb))
            *prod(factorial(2*c+1) for c in range(nu))
            *prod(factorial(2*nb+2*r+1) for r in range(nu)))
def sh1(n):
    nb=(n+1)//2
    return Lam0(n)*(n*(n+1)-1-4*nb)
def sh2(n):
    nb=(n+1)//2; nu=n//2
    return df(2*n-1)*df(2*nb-1)*2^nu*factorial(nu+1)*base(n)

_run("Corollary 27.15 (cor:03): Euler (s=0,t=3)",
     [("HH_n", lambda n: HH(n), base), ("HH_n^(1)", lambda n: HH(n,1), sh1), ("HH_n^(2)", lambda n: HH(n,2), sh2)])
