# check_Coro27_8.sage -- Corollary 27.8 (cor:gm1) of d51coro.tex: Gaussian c=-1.
# f=e^{-x+x^2/2}, a_n = (-1)^n I_n (signed involution numbers).
#   HH_n     = (-2)^C(n,2) prod_{k<n} k!
#   HH_n(f') = det(a_{2i+j+1}) = (-1)^n (-2)^C(n,2) prod_{k<n} k!
# Run:  sage check_Coro27_8.sage N
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

def df(k): return prod(range(k, 0, -2))
c = -1
def a(m): return sum(binomial(m, 2*j)*df(2*j-1)*c^(m-2*j) for j in range(m//2+1))
def HH(n, sh=0): return matrix(QQ, n, n, lambda i, j: a(2*i+j+sh)).det()
def RHS(n):  return (-2)^binomial(n,2)*prod(factorial(k) for k in range(1, n))
def RHS1(n): return (-1)^n*(-2)^binomial(n,2)*prod(factorial(k) for k in range(1, n))

_run("Corollary 27.8 (cor:gm1): Gaussian c=-1 (signed involution numbers)",
     [("HH_n", lambda n: HH(n), RHS), ("HH_n(f')", lambda n: HH(n,1), RHS1)])
