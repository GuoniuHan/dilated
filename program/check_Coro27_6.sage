# check_Coro27_6.sage -- Corollary 27.6 (cor:bqf) of d51coro.tex: (2m)!/m!.
# a_m = (2m)!/m!.  HH_n=det(a_{2i+j}), HH_n^(1)=det(a_{2i+j+1}).
#   HH_n   = 8^C(n,2) prod_i i!(4i)!/(2i)!
#   HH_n^1 = 8^C(n,2) prod_i i!(4i+2)!/(2i+1)!
# Run:  sage check_Coro27_6.sage N
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

def a(m): return factorial(2*m)/factorial(m)
def HH(n, sh=0): return matrix(QQ, n, n, lambda i, j: a(2*i+j+sh)).det()
def RHS(n):  return 8^binomial(n,2)*prod(factorial(i)*factorial(4*i)/factorial(2*i) for i in range(n))
def RHS1(n): return 8^binomial(n,2)*prod(factorial(i)*factorial(4*i+2)/factorial(2*i+1) for i in range(n))

_run("Corollary 27.6 (cor:bqf): a_m=(2m)!/m!",
     [("HH_n", lambda n: HH(n), RHS), ("HH_n^(1)", lambda n: HH(n,1), RHS1)])
