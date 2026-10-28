# Independent recomputation of the minimal obstructions for Phi(9).
# G = {C_n : n=1..34,36,39,40,42} u {C2 x C2n : n=1..13}; groups are pairs (m,n), m|n.
# Constraints: Weil pairing => phi(m) | 9 => m in {1,2}; primes dividing the order are <= 89.
from math import gcd
def primes_upto(N):
    return [p for p in range(2,N+1) if all(p%q for q in range(2,int(p**0.5)+1))]
def factor(n):
    f={};p=2
    while p*p<=n:
        while n%p==0: f[p]=f.get(p,0)+1; n//=p
        p+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
def phi(n):
    r=n
    for p in factor(n): r=r//p*(p-1)
    return r
# Weil pairing restriction
Ms=[m for m in range(1,200) if 9%phi(m)==0]
print("m with phi(m)|9 :",Ms)
cyc=set(list(range(1,35))+[36,39,40,42])
G=set((1,n) for n in cyc)|set((2,2*n) for n in range(1,14))
NVinf=set((1,n) for n in list(range(1,29))+[30,36])|set((2,2*n) for n in range(1,13))
print("#G =",len(G)," #Phi_inf(9) (NV Thm 1.1(3)) =",len(NVinf)," G\\Phi_inf =",sorted(G-NVinf))
def subgroups_types(m,n):
    # all (a,b) with a|b, a|m, b|n  (subgroup types of C_m x C_n with m|n)
    out=set()
    for a in range(1,m+1):
        if m%a: continue
        for b in range(1,n+1):
            if n%b==0 and b%a==0: out.add((a,b))
    return out
# closure check
for (m,n) in G:
    for s in subgroups_types(m,n):
        assert s in G, (m,n,s)
print("G closed under subgroups: True")
P=set(primes_upto(89))
def ok_primes(N): return all(p in P for p in factor(N))
mins=[]
# cyclic: bound n <= 42*89 suffices since n/p in cyc for all p|n; we search far beyond
for n in range(2,20000):
    if (1,n) in G or not ok_primes(n): continue
    if all((1,n//p) in G for p in factor(n)):
        # all proper subgroups
        assert all(s in G for s in subgroups_types(1,n) if s!=(1,n))
        mins.append((1,n))
for n in range(1,20000):
    if (2,2*n) in G or not ok_primes(2*n): continue
    if all(s in G for s in subgroups_types(2,2*n) if s!=(2,2*n)):
        mins.append((2,2*n))
cyc_m=[g for g in mins if g[0]==1]; nc_m=[g for g in mins if g[0]==2]
print("minimal obstructions: total",len(mins),"cyclic",len(cyc_m),"noncyclic",len(nc_m))
print("cyclic orders:",[g[1] for g in cyc_m])
print("noncyclic:",nc_m)
print("max cyclic order:",max(g[1] for g in cyc_m),"(bound 42*89 =",42*89,")")
# also: is every group outside G (m in {1,2}, primes<=89) containing a minimal one? by construction yes (finite poset)
# Breakdown of cyclic obstructions
pr=[n for (_,n) in cyc_m if len(factor(n))==1 and list(factor(n).values())[0]==1]
print("prime-order obstructions:",pr)
print("non-prime cyclic:",[n for (_,n) in cyc_m if n not in pr])
