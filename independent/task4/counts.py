# Closed-point counts on Y1(N) and Y1(2,2n) over F_p (p = 2,3,5) in degrees 1..9, computed from the weighted
# group-structure histograms produced by hist.gp (all F_q-isomorphism classes of elliptic curves, weight 1/#Aut).
#   #Y1(N)(F_q)   = sum_E #{P in E(F_q) of exact order N} / #Aut(E)          (N >= 4)
#   #Y1(2,2n)(F_q)= sum_E #{(P,Q): P of order 2, Q of order 2n, P not in <Q>} / #Aut(E)   (2n >= 4)
#   closed points of degree d: a_d = (1/d) sum_{e|d} mu(d/e) #Y(F_{p^e})
import glob, re, os, sys
from fractions import Fraction as Fr
from math import gcd
D = os.path.dirname(os.path.abspath(__file__))
H = {}   # (p,e) -> {(n1,n2): weight}
mass = {}
for fn in glob.glob(os.path.join(D, 'hist', 'h_*.txt')):
    txt = open(fn).read()
    assert 'CHUNK_DONE' in txt, fn
    for line in txt.splitlines():
        line = line.strip()
        if line.startswith('["MASS"'):
            v = line.strip('[]').split(','); p, e = int(v[1]), int(v[2]); mass[(p, e)] = mass.get((p, e), 0) + Fr(v[5].strip())
            continue
        if line.startswith('["CHUNK_DONE"'): continue
        v = [x.strip() for x in line.strip('[]').split(',')]
        p, e, k, n1, n2, w = int(v[0]), int(v[1]), int(v[2]), int(v[3]), int(v[4]), Fr(v[5])
        d = H.setdefault((p, e), {}); d[(n1, n2)] = d.get((n1, n2), 0) + w
for (p, e), m in sorted(mass.items()):
    assert m == p**e, (p, e, m)
def factor(n):
    f = {}; q = 2
    while q*q <= n:
        while n % q == 0: f[q] = f.get(q, 0) + 1; n //= q
        q += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f
def mu(n):
    f = factor(n)
    return 0 if any(v > 1 for v in f.values()) else (-1)**len(f)
def divisors(n): return [d for d in range(1, n+1) if n % d == 0]
def exact_order(N, n1, n2): return sum(mu(N//d) * gcd(d, n1) * gcd(d, n2) for d in divisors(N))
def npts(kind, N, p, e):
    tot = Fr(0)
    for (n1, n2), w in H[(p, e)].items():
        if kind == 'X1': c = exact_order(N, n1, n2)
        else: c = 2*exact_order(N, n1, n2) if n1 % 2 == 0 else 0   # N = 2n here
        tot += w * c
    assert tot.denominator == 1, (kind, N, p, e, tot)
    return int(tot)
def closed(kind, N, p, dmax=9):
    pts = {e: npts(kind, N, p, e) for e in range(1, dmax+1)}
    out = []
    for d in range(1, dmax+1):
        s = sum(mu(d//e) * pts[e] for e in divisors(d))
        assert s % d == 0
        out.append(s // d)
    return out
if __name__ == '__main__':
    # validation: X1(11) is the elliptic curve 11a3; #X1(11)(F_{p^e}) = p^e + 1 - (alpha^e + beta^e),
    # cusps: 5 rational + one degree-5 orbit (p = 2,3,5 have order 5 in (Z/11)^*/{+-1}).
    for p, ap in [(2, -2), (3, -1), (5, 1)]:
        s = [2, ap]
        for e in range(2, 10): s.append(ap*s[-1] - p*s[-2])
        for e in range(1, 10):
            full = p**e + 1 - s[e]
            cusps = 5 + (5 if e % 5 == 0 else 0)
            assert npts('X1', 11, p, e) == full - cusps, (p, e)
    print('VALIDATION X1(11) over F_p^e (p=2,3,5; e<=9) against 11a3 point counts: OK')
