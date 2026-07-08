# check_Coro27_44.sage -- Corollary 27.44 (cor:cxbeh) of d51coro.tex: mult-Bessel nu=1/2, f=(1+x)sin^2x/x^2.
#   HH_{2m} = (-1)^m 2^{m(8m-5)} prod_{k<m}(2k)!(2k+1)! prod_{i<2m} (2i)!i!(m-1+i)!/((2m+2i-1)!(m+i)!)
# Run:  sage check_Coro27_44.sage N   (checks HH_{2m} for m=1..N)
import sys
def poch(z,k): return prod(z+i for i in range(k))
nu = QQ(1)/2
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
def RHS(m): return ((-1)^m*2^(m*(8*m-5))*prod(factorial(2*k)*factorial(2*k+1) for k in range(m))
                    *prod(factorial(2*i)*factorial(i)*factorial(m-1+i)/(factorial(2*m+2*i-1)*factorial(m+i)) for i in range(2*m)))
print("Corollary 27.44 (cor:cxbeh): mult-Bessel nu=1/2, f=(1+x)sin^2x/x^2   [even orders]"); print()
ok=True
for m in range(1,N+1):
    d=HH2m(m)-RHS(m)
    if d!=0: ok=False
    print("n=2m=%d:  HH=%s (dif %s)"%(2*m, HH2m(m), d))
print(); print("Resume: all dif = 0. Formula checked." if ok else "Resume: some dif != 0. NOT checked!")
