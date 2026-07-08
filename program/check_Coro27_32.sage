# check_Coro27_32.sage -- Corollary 27.32 (elliptic m=2) of d51coro.tex.
# g_m = (1+sn(x,m))/cn(x,m);  HH_n(g_m) = (1-m)^C(n,2) HH_n(g_0),  HH_n(g_0)=2^{-C(n,2)} prod k!(2k)!.
# Run:  sage check_Coro27_32.sage N
import sys
def _run(title, claims, dN=6):
    N = int(sys.argv[1]) if len(sys.argv) > 1 else dN
    print(title); print(); ok = True
    for n in range(1, N+1):
        parts=[]
        for nm,L,R in claims:
            Lv=L(n); d=Lv-R(n)
            if d!=0: ok=False
            parts.append("%s=%s (dif %s)"%(nm,Lv,d))
        print("n=%d:  %s"%(n,"   ".join(parts)))
    print(); print("Resume: all dif = 0. Formula checked." if ok else "Resume: some dif != 0. NOT checked!")
xx=var('x')
def moments(Mmax):
    g=(1+jacobi('sn',xx,2))/jacobi('cn',xx,2)
    t=taylor(g,xx,0,Mmax)
    return [factorial(k)*QQ(t.coefficient(xx,k)) for k in range(Mmax+1)]
N0=int(sys.argv[1]) if len(sys.argv)>1 else 6
_A=moments(3*N0+2)
def a(m): return _A[m]
def HH(n): return matrix(QQ,n,n,lambda i,j: a(2*i+j)).det()
def base(n): return ((-1/2)^binomial(n,2))*prod(factorial(k)*factorial(2*k) for k in range(1,n))
_run("Corollary 27.32 (elliptic m=2): g_m=(1+sn)/cn",[("HH_n",lambda n:HH(n),base)])
