# check_Coro27_41.sage -- Corollary 27.41 (cor:cbessel01) of d51coro.tex: Bessel s=0,t=1/2.
# f = J_0(x) + Si(x) = cosb_0 + int cosb_{1/2}.  cosb_nu = sum_k (-1)^k x^{2k}/(4^k k! (nu+1)_k).
# a_{2k}=(2k)! [x^{2k}]cosb_s,  a_{2k+1}=(2k)! [x^{2k}]cosb_t.
#   HH_{2m}   = 2^{-C(2m,2)} prod_{k<m} ((2k)!)^2 (2m+2k)!/(4m+2k-1)!
#   HH_{2m+1} = (2m)!(4m)!/(2^{m(2m+5)} m!(3m)!) prod_{k<m} ((2k)!)^2 (2m+2k)!/(4m+2k+1)!
#   HH_n^(2)  = det(a_{2i+j+2}) = (-1)^n 3((2n)!)^2/(2^{n+1} n!(3n)!) HH_n
# Run:  sage check_Coro27_41.sage N
import sys
def _run(title, claims, dN=6):
    N = int(sys.argv[1]) if len(sys.argv) > 1 else dN
    print(title); print(); ok = True
    for n in range(1, N+1):
        parts=[]
        for nm,L,R in claims:
            Lv=L(n); d=Lv-R(n)
            if d!=0: ok=False
            parts.append("%s (dif %s)"%(nm,d))
        print("n=%d:  %s"%(n,"   ".join(parts)))
    print(); print("Resume: all dif = 0. Formula checked." if ok else "Resume: some dif != 0. NOT checked!")
def poch(z,k): return prod(z+i for i in range(k))
s, t = QQ(0), QQ(1)/2
def a(m):
    k = m//2; nu = s if m % 2 == 0 else t
    return factorial(2*k)*(-1)^k/(4^k*factorial(k)*poch(nu+1, k))
def HH(n, sh=0): return matrix(QQ, n, n, lambda i, j: a(2*i+j+sh)).det()
def base(n):
    if n % 2 == 0:
        m = n//2
        return 2^(-binomial(2*m,2))*prod(factorial(2*k)^2*factorial(2*m+2*k)/factorial(4*m+2*k-1) for k in range(m))
    m = (n-1)//2
    return (factorial(2*m)*factorial(4*m)/(2^(m*(2*m+5))*factorial(m)*factorial(3*m))
            *prod(factorial(2*k)^2*factorial(2*m+2*k)/factorial(4*m+2*k+1) for k in range(m)))
def sh2(n): return (-1)^n*3*factorial(2*n)^2/(2^(n+1)*factorial(n)*factorial(3*n))*base(n)
_run("Corollary 27.41 (cor:cbessel01): Bessel s=0,t=1/2, f=J_0+Si",
     [("HH_n", lambda n: HH(n), base), ("HH_n^(2)", lambda n: HH(n,2), sh2)])
