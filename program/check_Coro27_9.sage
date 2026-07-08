# check_Coro27_9.sage -- Corollary 27.9 (cor:euler) of d51coro.tex: Euler numbers (s=0,t=1).
# a_{2k}=E^{(s+1)}_{2k} (S-fraction u_j=j(j+s)), a_{2k+1}=E^{(t+1)}_{2k} (v_j=j(j+t)); s=0,t=1.
#   HH_n     = 2^{-C(n,2)} prod_{k<n} k!(2k)!
#   HH_n^(1) = det(a_{2i+j+1}) = (2n-1)!! HH_n
#   HH_n^(2) = det(a_{2i+j+2}) = (2n-1)!!(2nbar-1)!! 2^nund nund! HH_n
# Run:  sage check_Coro27_9.sage N
import sys
def _run(title, claims, dN=6):
    N = int(sys.argv[1]) if len(sys.argv) > 1 else dN
    print(title); print(); ok = True
    for n in range(1, N+1):
        parts = []
        for nm, L, R in claims:
            Lv = L(n); d = Lv - R(n)
            if d != 0: ok = False
            parts.append("%s (dif %s)" % (nm, d))
        print("n=%d:  %s" % (n, "   ".join(parts)))
    print(); print("Resume: all dif = 0. Formula checked." if ok else "Resume: some dif != 0. NOT checked!")

s, t = 0, 1
def df(k): return prod(range(k, 0, -2))
def sfrac(bfun, K):
    PS = PowerSeriesRing(QQ, 'x', default_prec=K+2); x = PS.gen(); f = PS(1)
    for j in range(K+1, 0, -1): f = (1 - bfun(j)*x*f)^(-1)
    return [f[k] for k in range(K+1)]
_MU = None
def a(m):
    global _MU
    return _MU[0][m//2] if m % 2 == 0 else _MU[1][m//2]
def HH(n, sh=0): return matrix(QQ, n, n, lambda i, j: a(2*i+j+sh)).det()
def setup(N):
    global _MU
    K = 2*N + 2
    _MU = (sfrac(lambda j: j*(j+s), K), sfrac(lambda j: j*(j+t), K))

def base(n): return 2^(-binomial(n,2))*prod(factorial(k)*factorial(2*k) for k in range(1, n))
def sh1(n):  return df(2*n-1)*base(n)
def sh2(n):
    nb=(n+1)//2; nu=n//2
    return df(2*n-1)*df(2*nb-1)*2^nu*factorial(nu)*base(n)

setup(int(sys.argv[1]) if len(sys.argv) > 1 else 6)
_run("Corollary 27.9 (cor:euler): Euler numbers (s=0,t=1)",
     [("HH_n", lambda n: HH(n), base),
      ("HH_n^(1)", lambda n: HH(n,1), sh1),
      ("HH_n^(2)", lambda n: HH(n,2), sh2)])
