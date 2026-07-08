# check_Coro27_43.sage -- Corollary 27.43 (cor:cxbe0) of d51coro.tex: mult-Bessel nu=0, f=(1+x)J_0^2.
# cosb_nu = sum_k (-1)^k x^{2k}/(4^k k!(nu+1)_k); f=(1+x)cosb_nu^2; a_n=n![x^n]f.  Even orders only.
#   HH_{2m} = (-1)^m/2^{m(6m-5)} prod_{k<m}((2k)!/k!)^4 prod_{i<2m} ((2i)!)^2/(i!((m-1+i)!)^2)
# Run:  sage check_Coro27_43.sage N   (checks HH_{2m} for m=1..N)
import sys
def poch(z,k): return prod(z+i for i in range(k))
nu = QQ(0)
def moments(Mmax):
    prec=Mmax+1
    PS=PowerSeriesRing(QQ,'x',default_prec=prec+1); x=PS.gen()
    cb=sum((-1)^k*x^(2*k)/(4^k*factorial(k)*poch(nu+1,k)) for k in range(prec//2+1))
    f=(1+x)*cb^2
    return [factorial(m)*f[m] for m in range(Mmax+1)]
N=int(sys.argv[1]) if len(sys.argv)>1 else 4
_A=moments(6*N+2)
def a(m): return _A[m]
def HH2m(m): return matrix(QQ,2*m,2*m,lambda i,j: a(2*i+j)).det()
def RHS(m): return ((-1)^m/2^(m*(6*m-5))*prod((factorial(2*k)/factorial(k))^4 for k in range(m))
                    *prod(factorial(2*i)^2/(factorial(i)*factorial(m-1+i)^2) for i in range(2*m)))
print("Corollary 27.43 (cor:cxbe0): mult-Bessel nu=0, f=(1+x)J_0^2   [even orders]"); print()
ok=True
for m in range(1,N+1):
    d=HH2m(m)-RHS(m)
    if d!=0: ok=False
    print("n=2m=%d:  HH=%s (dif %s)"%(2*m, HH2m(m), d))
print(); print("Resume: all dif = 0. Formula checked." if ok else "Resume: some dif != 0. NOT checked!")
