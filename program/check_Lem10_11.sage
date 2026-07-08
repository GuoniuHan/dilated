# check_Lem10_11.sage -- verification of Lemma 10.11 (lem:t3contig) of d05eulershift.tex.
#
# Contiguous relation in s at t=3.
#   (eq:t3moment)  a_{2k}(s) = [ a_{2k+2}(s-2) + (s-1)^2 a_{2k}(s-2) ] / (s(s-1)),
#                  a_{2k+1}(s) = a_{2k+1}(s-2)   (odd part is s-independent).
#   (eq:t3detcontig)  (s(s-1))^{nund} HH_n^{(1)}(s) = sum_{q=0}^{nund} (s-1)^{2q} G_q(s-2),
#   where G_q is the minor at parameter s-2 whose even columns are unchanged and
#   whose nund odd columns carry the even-moment bases {1,...,nund+1}\{q+1}.
#
# s symbolic (t=3).  Part A checks the moment relation; Part B the determinant one.
#
# Run:  sage check_Lem10_11.sage N

import sys

R = PolynomialRing(QQ, 's'); s = R.gen(); F = R.fraction_field()

def cS(i):   return 2*i*(2*i+s) + (2*i+1)*(2*i+1+s)
def lamS(i): return (2*i-1)*(2*i-1+s) * 2*i*(2*i+s)
def cT(m):   return 2*m*(2*m+3) + (2*m+1)*(2*m+4)
def lamT(m): return (2*m-1)*(2*m+2) * 2*m*(2*m+3)

def moments(cf, lf, Kmax):
    sz = Kmax + 1
    J = matrix(F, sz, sz)
    for i in range(sz):
        J[i,i] = cf(i)
        if i+1 < sz: J[i,i+1] = F(1); J[i+1,i] = lf(i+1)
    mu = [F(1)]; Jp = identity_matrix(F, sz)
    for k in range(1, Kmax+1): Jp = Jp*J; mu.append(Jp[0,0])
    return mu

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 5

    K = 3*N + 3
    muS = moments(cS, lamS, K)   # a_{2k} = muS[k]  (depends on s)
    muT = moments(cT, lamT, K)   # a_{2k+1} = muT[k] (s-independent)

    def sm2(expr): return F(expr).subs(s=s-2)   # value at parameter s-2

    all_zero = True

    # ---- Part A: eq:t3moment ----
    print("Lemma 10.11 (lem:t3contig) of d05eulershift.tex   [s symbolic, t=3]")
    print("")
    print("Part A -- eq:t3moment (moment recurrence in s)")
    print("  k | dif even = a_{2k}(s) - [a_{2k+2}(s-2)+(s-1)^2 a_{2k}(s-2)]/(s(s-1)) | dif odd")
    print("----+---------------------------------------------------------------------+--------")
    for k in range(2*N + 1):
        de = F(muS[k] - (sm2(muS[k+1]) + (s-1)^2*sm2(muS[k])) / (s*(s-1)))
        do = F(muT[k] - sm2(muT[k]))
        if de != 0 or do != 0: all_zero = False
        print("%3d | %67s | %s" % (k, de, do))

    # ---- Part B: eq:t3detcontig ----
    def a(m, at_s):
        # a_m at parameter 'at_s'; even m uses muS (subs), odd m uses muT (const)
        return (F(muS[m//2]).subs(s=at_s) if m % 2 == 0 else muT[m//2])

    print("")
    print("Part B -- eq:t3detcontig ((s(s-1))^{nund} HH_n^{(1)}(s) = sum_q (s-1)^{2q} G_q(s-2))")
    print("  n | dif = LHS - RHS")
    print("----+----------------")
    for n in range(1, N + 1):
        nund = n//2
        HH = matrix(F, n, n, lambda i, j: a(2*i + j + 1, s)).det()
        LHS = F((s*(s-1))^nund * HH)

        RHS = F(0)
        for q in range(nund + 1):
            bases = [b for b in range(1, nund + 2) if b != q + 1]  # {1..nund+1}\{q+1}
            def entry(i, j):
                if j % 2 == 0:                       # even column: odd moment, s-indep
                    return a(2*i + j + 1, s - 2)
                base = bases[(j - 1)//2]             # odd column -> even-moment base
                return F(muS[i + base]).subs(s=s-2)
            Gq = matrix(F, n, n, entry).det()
            RHS += (s-1)^(2*q) * Gq
        d = F(LHS - RHS)
        if d != 0: all_zero = False
        print("%3d | %s" % (n, d))

    print("")
    if all_zero:
        print("Resume: all dif = 0. Lemma checked.")
    else:
        print("Resume: some dif != 0. Lemma NOT checked!")
