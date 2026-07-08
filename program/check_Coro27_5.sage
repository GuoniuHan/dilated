# check_Coro27_5.sage -- Corollary 27.5 (cor:bddf) of d51coro.tex: double factorials.
# Two base sequences: a1_m=(2m-1)!!, a2_m=(2m+1)!!.  df(k)=k!! (empty product=1).
#   HH_n((2m-1)!!)   = 4^C(n,2) prod_i i!(4i-1)!!
#   HH_n((2m+1)!!)   = 4^C(n,2) prod_i i!(4i+1)!!
#   HH_n^(1)((2m+1)!!)= 4^C(n,2) prod_i i!(4i+3)!!
# Run:  sage check_Coro27_5.sage N
import sys
def _run(title, claims, dN=6):
    N = int(sys.argv[1]) if len(sys.argv) > 1 else dN
    print(title); print(); ok = True
    for n in range(1, N+1):
        parts = []
        for nm, L, R in claims:
            Lv = L(n); d = Lv - R(n)
            if d != 0: ok = False
            parts.append("%s (dif %s)" % (nm, d))
        print("n=%d:  %s" % (n, "   ".join(parts)))
    print(); print("Resume: all dif = 0. Formula checked." if ok else "Resume: some dif != 0. NOT checked!")

def df(k): return prod(range(k, 0, -2))          # k!!, empty product = 1 (so (-1)!!=1)
def a1(m): return df(2*m-1)
def a2(m): return df(2*m+1)
def HH(a, n, sh=0): return matrix(QQ, n, n, lambda i, j: a(2*i+j+sh)).det()

_run("Corollary 27.5 (cor:bddf): double factorials",
     [("HH((2m-1)!!)",  lambda n: HH(a1,n),  lambda n: 4^binomial(n,2)*prod(factorial(i)*df(4*i-1) for i in range(n))),
      ("HH((2m+1)!!)",  lambda n: HH(a2,n),  lambda n: 4^binomial(n,2)*prod(factorial(i)*df(4*i+1) for i in range(n))),
      ("HH1((2m+1)!!)", lambda n: HH(a2,n,1),lambda n: 4^binomial(n,2)*prod(factorial(i)*df(4*i+3) for i in range(n)))])
