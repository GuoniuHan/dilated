# check_Coro27_1.sage -- Corollary 27.1 (cor:bcatalan) of d51coro.tex: Catalan numbers.
# a_m = C_m = binomial(2m,m)/(m+1).  HH_n=det(a_{2i+j}), HH_n^(1)=det(a_{2i+j+1}).
#   HH_n   = 2^C(n,2) prod_i binom(4i,2i) i!(n+i)!/(n+2i)!             = 1,3,32,1232,...
#   HH_n^1 = 2^C(n,2) prod_i i!(4i+2)!(n+i)!/((2i)!(2i+1)!(n+2i+1)!)
# Run:  sage check_Coro27_1.sage N
import sys
def _run(title, claims, dN=6):
    N = int(sys.argv[1]) if len(sys.argv) > 1 else dN
    print(title); print()
    ok = True
    for n in range(1, N+1):
        parts = []
        for nm, L, R in claims:
            Lv = L(n); d = Lv - R(n)
            if d != 0: ok = False
            parts.append("%s=%s (dif %s)" % (nm, Lv, d))
        print("n=%d:  %s" % (n, "   ".join(parts)))
    print()
    print("Resume: all dif = 0. Formula checked." if ok else "Resume: some dif != 0. NOT checked!")

def a(m): return binomial(2*m, m)/(m+1)
def HH(n, sh=0): return matrix(QQ, n, n, lambda i, j: a(2*i+j+sh)).det()
def RHS(n):  return 2^binomial(n,2)*prod(binomial(4*i,2*i)*factorial(i)*factorial(n+i)/factorial(n+2*i) for i in range(n))
def RHS1(n): return 2^binomial(n,2)*prod(factorial(i)*factorial(4*i+2)*factorial(n+i)/(factorial(2*i)*factorial(2*i+1)*factorial(n+2*i+1)) for i in range(n))

_run("Corollary 27.1 (cor:bcatalan): Catalan numbers  a_m=C_m",
     [("HH_n", lambda n: HH(n), RHS), ("HH_n^(1)", lambda n: HH(n,1), RHS1)])
