# r1root.sage -- self-contained
#
# Roots of the Poupard and Kreweras polynomials (Chapoton-Han conjecture)
#
# Builds the rho-matrix M whose entry (i,m) is
#     rho( x^i (x+1)^m ) / 2^m * 2^i ,
# (even and odd powers m interleaved), and its row-difference matrix A.
#
# rho is the "constant-term" reduction: repeatedly apply the operator ND
# (ParamD=0) until the polynomial has degree <= 1, then read off its value
# at x=0.
# 

def ND(P, ParamD=2):
    d = P.degree(x) + P.low_degree(x)
    T = (x^(d + 2*ParamD) + 1)*P(x=1) - 2*x^ParamD*P
    return expand(factor(T/(x-1)^2))


def rho(f):
    P = f
    while P.degree(x) + P.low_degree(x) >= 2:
        P = ND(P, ParamD=0)
    return P(x=0)


# ---- matrix M ----
M = []
sz = 6

L1 = []
i = 0
for m in range(3*sz/2):
    m2 = 2*m
    L1.append(rho(x^i*(x+1)^m2)/2^m*2^i)
    m2 = 2*m + 1
    L1.append(rho(x^i*(x+1)^m2)/2^m*2^i)
print("The first row:")
print(L1)

for i in range(sz):
    L = []
    for m in range(sz/2):
        m2 = 2*m
        L.append(rho(x^i*(x+1)^m2)/2^m*2^i)
        m2 = 2*m + 1
        L.append(rho(x^i*(x+1)^m2)/2^m*2^i)
    M.append(L)

print()
print("The modified Chapoton-Han matrix (the k-th row multiplied by 2^k):")
print(matrix(M))
#print(factor(Integer(det(matrix(M)))))


# ---- row-difference matrix A: line(i) - line(i-1) ----
print()
print("First difference: row(k) - row(k-1) for k = 1, 2, 3, ...:")
sz = len(M)
A = [copy(m) for m in M]
for i in range(1, sz):
    for j in range(sz):
        A[i][j] = M[i][j] - M[i-1][j]
print(matrix(A))
#print(factor(Integer(det(matrix(A)))))

print()
print("Second difference: row(k) - row(k-1) for k = 2, 3, ...:")
print("...")
print("This yields the dilated matrix.")




