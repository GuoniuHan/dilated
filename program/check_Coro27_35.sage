# check_Coro27_35.sage -- Corollary 27.35 (algebraic s=-1) of d51coro.tex.
# f=(1+x)/(1-x^2)^{s/2}: a_{2k}=(2k)!/k! (s/2)_k, a_{2k+1}=(2k+1)a_{2k}.
#   HH_n = (prod_{k<n}(2k)!) prod_{j=1}^{n-1}(s+2j-2)^{n-j}.
# Run:  sage check_Coro27_35.sage N
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
def poch(z,k): return prod(z+i for i in range(k))
sval=-1
def a(m):
    k=m//2
    a2=factorial(2*k)/factorial(k)*poch(QQ(sval)/2,k)
    return a2 if m%2==0 else (2*k+1)*a2
def HH(n): return matrix(QQ,n,n,lambda i,j: a(2*i+j)).det()
def base(n): return prod(factorial(2*k) for k in range(1,n))*prod((sval+2*j-2)^(n-j) for j in range(1,n))
_run("Corollary 27.35 (algebraic s=-1): f=(1+x)/(1-x^2)^{s/2}",[("HH_n",lambda n:HH(n),base)])
