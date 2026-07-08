# check_Coro27_2.sage -- Corollary 27.2 (cor:bcbinom) of d51coro.tex: central binomial coeffs.
# a_m = binomial(2m,m).  HH_n=det(a_{2i+j}), HH_n^(1)=det(a_{2i+j+1}).
#   HH_n   = 2^{C(n,2)+n-1} prod_i binom(4i,2i) i!(n-1+i)!/(n-1+2i)!        = 1,8,224,...
#   HH_n^1 = 2^{C(n,2)+n-1} prod_i i!(4i+2)!(n-1+i)!/((2i)!(2i+1)!(n+2i)!)
# Run:  sage check_Coro27_2.sage N
import sys
def _run(title, claims, dN=6):
    N = int(sys.argv[1]) if len(sys.argv) > 1 else dN
    print(title); print(); ok = True
    for n in range(1, N+1):
        parts = []
        for nm, L, R in claims:
            Lv = L(n); d = Lv - R(n)
            if d != 0: ok = False
            parts.append("%s=%s (dif %s)" % (nm, Lv, d))
        print("n=%d:  %s" % (n, "   ".join(parts)))
    print(); print("Resume: all dif = 0. Formula checked." if ok else "Resume: some dif != 0. NOT checked!")

def a(m): return binomial(2*m, m)
def HH(n, sh=0): return matrix(QQ, n, n, lambda i, j: a(2*i+j+sh)).det()
def RHS(n):  return 2^(binomial(n,2)+n-1)*prod(binomial(4*i,2*i)*factorial(i)*factorial(n-1+i)/factorial(n-1+2*i) for i in range(n))
def RHS1(n): return 2^(binomial(n,2)+n-1)*prod(factorial(i)*factorial(4*i+2)*factorial(n-1+i)/(factorial(2*i)*factorial(2*i+1)*factorial(n+2*i)) for i in range(n))

_run("Corollary 27.2 (cor:bcbinom): central binomial coefficients  a_m=binom(2m,m)",
     [("HH_n", lambda n: HH(n), RHS), ("HH_n^(1)", lambda n: HH(n,1), RHS1)])
