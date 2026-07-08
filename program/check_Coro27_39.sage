# check_Coro27_39.sage -- Corollary 27.39 (cor:calgsq2) of d51coro.tex: sq-algebraic s=2.
# f=(1+x)^2/(1-x^2)=(1+x)/(1-x): a_0=1, a_m=2 m! (m>=1).
#   HH_n = 2^C(n,2) prod_{k<n}(2k)! prod_{i=1}^{n-2} i! (2^n (n-1)! - (2n-1)!!) = 1,4,384,-39813120,...
# Run:  sage check_Coro27_39.sage N
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
def df(k): return prod(range(k, 0, -2))
def a(m): return QQ(1) if m==0 else 2*factorial(m)
def HH(n): return matrix(QQ,n,n,lambda i,j: a(2*i+j)).det()
def base(n): return 2^binomial(n,2)*prod(factorial(2*k) for k in range(1,n))*prod(factorial(i) for i in range(1,n-1))*(2^n*factorial(n-1)-df(2*n-1))
_run("Corollary 27.39 (cor:calgsq2): sq-algebraic s=2, f=(1+x)/(1-x)",[("HH_n",lambda n:HH(n),base)])
