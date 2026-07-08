## =====================================================================
## r2coro.sage
## Reproduces every Corollary of d51coro.tex, in the same order, printing
## the dilated Hankel determinants HH_n = det(a_{2i+j}) (and, where the
## corollary records them, the single / double shifts sHH / ssHH).
##
## Machinery (hankel, dhankel, hfactor, int_sec, the taylor->sequence
## extraction) is copied from github/s10famst.sage and s11xcosx.sage.
## The Bessel cosine cosb_nu = Gamma(nu+1)(2/x)^nu J_nu(x) is built as a
## truncated power series, as in program/d24proof.sage.
## Output format is the same as in those programs.
## =====================================================================

def hfactor(a):
    if a==0: return a
    return factor(a)

def hankel(s, m): # <<<
  A=[]
  L=s.copy()
  A.append(L)
  for j in range(m-1):
    U=[L[i+1] for i in range(len(L)-1)]
    A.append(U)
    L=U.copy()
  B=[]
  for a in A:
    B.append(a[:m])
  return(matrix(B))
# >>>

def dhankel(s, m): # <<<
  A=[]
  L=s.copy()
  A.append(L)
  for j in range(m-1):
    U=[L[i+2] for i in range(len(L)-2)]
    A.append(U)
    L=U.copy()
  B=[]
  for a in A:
    B.append(a[:m])
  return(matrix(B))
# >>>

def int_sec(n, x): # <<<  integral of sec(x)^n dx  (copied from s10)
    if n < 0:
        m=-n
        return integral(cos(x)^m, x)
    if n == 0:
        return x
    if n == 1:
        return log(sec(x) + tan(x))
    return (sec(x)^(n-2)*tan(x))/(n-1) \
         + (n-2)/(n-1)*int_sec(n-2, x)
# >>>

## -- global size of the determinant table ---------------------------
M = 8

def egf_seq(f, N): # exponential generating function -> moment sequence
    p = taylor(f, x, 0, N+2)
    return [QQ(p.coefficient(x,j)*factorial(j)) for j in range(N+1)]

def print_dets(seq, tag, m):
    for j in range(1, m+1):
        da = hfactor(det(dhankel(seq, j)))
        print("%s%d ="%(tag, j), da)

def report(label, f, shifts=(), m=M):
    """label + egf f; shifts subset of {1,2} for single/double shift."""
    Deg = 3*m + 3
    print()
    print(label)
    print("case:coro| f=", f)
    seq = egf_seq(f, Deg)
    print("case:coro| seq=", seq[:14])
    print_dets(seq, "HH", m)
    if 1 in shifts:
        f1 = diff(f, x)
        print("case:coro| f'=", f1)
        print_dets(egf_seq(f1, Deg), "sHH", m)
    if 2 in shifts:
        f2 = diff(f, x, 2)
        print("case:coro| f''=", f2)
        print_dets(egf_seq(f2, Deg), "ssHH", m)

def betaseq(r, a, b, N): # plain sequence a_m = r^m (a)_m / (b)_m
    s = [QQ(1)]
    for mm in range(N):
        ratio = r*(mm + a)
        if b is not None:
            ratio = ratio/(mm + b)
        s.append(s[-1]*ratio)
    return [QQ(v) for v in s]

def report_seq(label, fdesc, seq, shift=False, m=M):
    print()
    print(label)
    print("case:coro| f=", fdesc)
    print("case:coro| seq=", seq[:14])
    print_dets(seq, "HH", m)
    if shift:
        print_dets(seq[1:], "sHH", m)

## -- Bessel cosine cosb_nu and its moment sequences -----------------
def poch(a, k):
    return prod(a + j for j in range(k))

def bessel_moments(s, t, N):    # egf moments of  cosb_s + int_0^x cosb_t
    P = PowerSeriesRing(QQ, 'x', default_prec=N+3); x = P.gen()
    cs = sum((-1)^k/(4^k*factorial(k)*poch(s+1,k))*x^(2*k)
             for k in range((N+2)//2+1))
    ct = sum((-1)^k/(4^k*factorial(k)*poch(t+1,k))/(2*k+1)*x^(2*k+1)
             for k in range((N+2)//2+1))
    f = cs + ct
    return [QQ(factorial(k)*f[k]) for k in range(N+1)]

def xbessel_moments(nu, N):     # egf moments of  (1+x) cosb_nu^2
    P = PowerSeriesRing(QQ, 'x', default_prec=N+3); x = P.gen()
    cb = sum((-1)^k/(4^k*factorial(k)*poch(nu+1,k))*x^(2*k)
             for k in range((N+2)//2+1))
    f = (1 + x)*cb^2
    return [QQ(factorial(k)*f[k]) for k in range(N+1)]

def report_bessel(label, fdesc, s, t, dshift=False, m=M):
    Deg = 3*m + 3
    print()
    print(label)
    print("case:coro| f=", fdesc)
    seq = bessel_moments(s, t, Deg)
    print("case:coro| seq=", seq[:10])
    print_dets(seq, "HH", m)
    if dshift:                  # double shift = a_{k+2}
        print_dets(seq[2:], "ssHH", m)

def report_xbessel(label, fdesc, nu, m=M):
    Deg = 3*m + 3
    print()
    print(label)
    print("case:coro| f=", fdesc)
    seq = xbessel_moments(nu, Deg)
    print("case:coro| seq=", seq[:10])
    print_dets(seq, "HH", m) # only even orders admit a product form


## =====================================================================
## 1.  The Beta family   a_m = r^m (a)_m / (b)_m   (with single shift)
## =====================================================================
print("#### The Beta family a_m = r^m (a)_m/(b)_m ####")
Nb = 3*M + 4
# Catalan  (4,1/2,2)
report_seq("cor:bcatalan | [Beta] Catalan (r,a,b)=(4,1/2,2)",
           "a_m = 4^m (1/2)_m/(2)_m = C(2m,m)/(m+1)  (Catalan)",
           betaseq(4, QQ(1)/2, 2, Nb), shift=True)
# central binomial  (4,1/2,1)
report_seq("cor:bcbinom | [Beta] central binomial (4,1/2,1)",
           "a_m = 4^m (1/2)_m/(1)_m = C(2m,m)",
           betaseq(4, QQ(1)/2, 1, Nb), shift=True)
# binomial(2m+1,m)  (4,3/2,2)
report_seq("cor:b2m1 | [Beta] binomial(2m+1,m) (4,3/2,2)",
           "a_m = 4^m (3/2)_m/(2)_m = C(2m+1,m)",
           betaseq(4, QQ(3)/2, 2, Nb), shift=True)
# factorials  (1,1,-)
report_seq("cor:bfact | [Beta] factorials m! (1,1,-)",
           "a_m = (1)_m = m!",
           betaseq(1, 1, None, Nb), shift=True)
# double factorials: (2m-1)!! (2,1/2,-),  (2m+1)!! (2,3/2,-) with its shift
report_seq("cor:bddf | [Beta] double factorial (2m-1)!! (2,1/2,-)",
           "a_m = 2^m (1/2)_m = (2m-1)!!",
           betaseq(2, QQ(1)/2, None, Nb), shift=False)
report_seq("cor:bddf | [Beta] double factorial (2m+1)!! (2,3/2,-)  (+ shift)",
           "a_m = 2^m (3/2)_m = (2m+1)!!",
           betaseq(2, QQ(3)/2, None, Nb), shift=True)
# (2m)!/m!  (4,1/2,-)
report_seq("cor:bqf | [Beta] (2m)!/m! (4,1/2,-)",
           "a_m = 4^m (1/2)_m = (2m)!/m!",
           betaseq(4, QQ(1)/2, None, Nb), shift=True)


## =====================================================================
## 2.  The Gaussian family            f = e^{c x + x^2/2}   (single shift)
## =====================================================================
print()
print("#### The Gaussian family ####")
for (lab, c) in [("cor:g1", 1), ("cor:gm1", -1)]:
    f = exp(c*x + x^2/2)
    report("%s | [Gauss]  c = %s"%(lab, c), f, shifts=(1,))


## =====================================================================
## 3.  The Euler number family        f = 1/cos^{s+1} + int_sec(t+1)
##     order and recorded shifts exactly as in d51coro.tex
## =====================================================================
print()
print("#### The Euler number family ####")
# (label, s, t, shifts):  single shift only for s in {0,2}, t in {1,3}
euler_cases = [
    ("cor:euler", 0, 1, (1, 2)),
    ("cor:12",    1, 2, (2,)),
    ("cor:23",    2, 3, (1, 2)),
    ("cor:10",    1, 0, (2,)),
    ("cor:21",    2, 1, (1, 2)),
    ("cor:32",    3, 2, (2,)),
    ("cor:03",    0, 3, (1, 2)),
    ("cor:30",    3, 0, (2,)),
]
for (lab, s, t, sh) in euler_cases:
    f = 1/cos(x)^(s+1) + int_sec(t+1, x)
    report("%s | [Euler]  [s, t] = %s"%(lab, [s, t]), f, shifts=sh)


## =====================================================================
## 4.  The secant-number family       f = (1+x)/cos^{s+1}
## =====================================================================
print()
print("#### The secant-number family (1+x)/cos^{s+1} ####")
xc_cases = [("cor:xc0", 0, (1,)), ("cor:xc1", 1, (1, 2)),
            ("cor:xc2", 2, (1,)), ("cor:xc3", 3, (1,))]
for (lab, s, sh) in xc_cases:
    f = (1 + x)/cos(x)^(s+1)
    report("%s | [xcos]  s = %s"%(lab, s), f, shifts=sh)


## =====================================================================
## 5.  Rank-one perturbation of secant-tangent
##     f_s = (sin x + 1)/cos^2 x + s sin x     (d10 = s=-1)
## =====================================================================
print()
print("#### Rank-one perturbation (sin x + 1)/cos^2 x + s sin x ####")
sin_cases = [("cor:sinm2", -2), ("cor:sin3", -1), ("cor:sin0", 0),
             ("cor:sin1", 1), ("cor:sin2p", 2)]
for (lab, s) in sin_cases:
    f = (sin(x) + 1)/cos(x)^2 + s*sin(x)
    report("%s | [sin_s]  s = %s"%(lab, s), f)


## =====================================================================
## 6.  The Springer secant ladder   f = 1/(cos x - t sin x)^r,  t=1
## =====================================================================
print()
print("#### The Springer secant ladder 1/(cos x - sin x)^r ####")
for (lab, r) in [("cor:spr1", 1), ("cor:spr2", 2)]:
    f = 1/(cos(x) - sin(x))^r
    report("%s | [Springer]  t=1, r = %s"%(lab, r), f)


## =====================================================================
## 7.  The derivative ladder   f = (cos x + sin x)/(cos x - sin x)^s
## =====================================================================
print()
print("#### The derivative ladder (cos x + sin x)/(cos x - sin x)^s ####")
for (lab, s) in [("cor:der1", 1), ("cor:der2", 2), ("cor:der3", 3)]:
    f = (cos(x) + sin(x))/(cos(x) - sin(x))^s
    report("%s | [Deriv]  s = %s"%(lab, s), f)


## =====================================================================
## 8.  The reciprocal-sine case   f = (1+x) x / sin x   (no free parameter)
## =====================================================================
print()
print("#### The reciprocal-sine case (1+x) x / sin x ####")
report("cor:xsin | [xsin]  (1+x) x / sin x", (1 + x)*x/sin(x))


## =====================================================================
## 9.  Elliptic deformation of the Euler numbers   g = (1+sn)/cn
## =====================================================================
print()
print("#### The elliptic deformation (1+sn(x,m))/cn(x,m) ####")
for (lab, mod) in [("cor:ellm1", -1), ("cor:ell2", 2), ("cor:ellh", QQ(1)/2)]:
    g = (1 + jacobi('sn', x, mod))/jacobi('cn', x, mod)
    report("%s | [Elliptic]  m = %s"%(lab, mod), g)


## =====================================================================
## 10.  The algebraic family          f = (1+x)/(1-x^2)^{s/2}
## =====================================================================
print()
print("#### The algebraic family (1+x)/(1-x^2)^{s/2} ####")
alg_cases = [("cor:algm3", -3), ("cor:algm1", -1), ("cor:alg1", 1),
             ("cor:alg2", 2), ("cor:alg3", 3)]
for (lab, s) in alg_cases:
    f = (1 + x)*(1 - x^2)^(-QQ(s)/2)
    report("%s | [alg]  s = %s"%(lab, s), f)


## =====================================================================
## 11.  The squared algebraic family   f = (1+x)^2/(1-x^2)^{s/2},  s != 3
## =====================================================================
print()
print("#### The squared algebraic family (1+x)^2/(1-x^2)^{s/2} ####")
for (lab, s) in [("cor:calgsq2", 2), ("cor:calgsq4", 4)]:
    f = (1 + x)^2*(1 - x^2)^(-QQ(s)/2)
    report("%s | [algsq]  s = %s"%(lab, s), f)


## =====================================================================
## 12.  The Bessel (s,t) family   f = cosb_s + int_0^x cosb_t   (double shift)
## =====================================================================
print()
print("#### The Bessel (s,t) family cosb_s + int_0^x cosb_t ####")
report_bessel("cor:cbessel01 | [Bessel]  (s,t) = (0, 1/2)  = J_0 + Si",
              "cosb_0 + int cosb_{1/2} = J_0(x) + Si(x)",
              QQ(0), QQ(1)/2, dshift=True)
report_bessel("cor:cbessel10 | [Bessel]  (s,t) = (1/2, 0)  = sin x/x + int J_0",
              "cosb_{1/2} + int cosb_0 = sin(x)/x + int_0^x J_0",
              QQ(1)/2, QQ(0), dshift=True)


## =====================================================================
## 13.  The multiplicative Bessel family   f = (1+x) cosb_nu^2  (even orders)
## =====================================================================
print()
print("#### The multiplicative Bessel family (1+x) cosb_nu^2 ####")
report_xbessel("cor:cxbe0 | [xBessel]  nu = 0    = (1+x) J_0^2",
               "(1+x) cosb_0^2 = (1+x) J_0(x)^2", QQ(0))
report_xbessel("cor:cxbeh | [xBessel]  nu = 1/2  = (1+x) sin^2 x / x^2",
               "(1+x) cosb_{1/2}^2 = (1+x) sin^2(x)/x^2", QQ(1)/2)
report_xbessel("cor:cxbe1 | [xBessel]  nu = 1    = 4(1+x) J_1^2 / x^2",
               "(1+x) cosb_1^2 = 4(1+x) J_1(x)^2/x^2", QQ(1))
