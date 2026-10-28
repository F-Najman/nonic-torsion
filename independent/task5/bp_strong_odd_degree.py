# Independent derivation of G_CM(9) from Bourdon-Pollack, "Torsion subgroups of CM elliptic curves over
# odd degree number fields" (arXiv:1601.00351v1), Theorem 1.1 (Odd Degree Theorem) + Theorem 1.2.
# Class numbers h(-l) computed with PARI (qfbclassno) for primes l = 3 mod 4.
import cypari2
pari = cypari2.Pari()
def delta(l, n):
    if l > 3: return (3*n)//2 - 1
    return 0 if n == 1 else (3*n)//2 - 2
def groups(d, lmax=2000, nmax=8):
    G = ['C1', 'C2', 'C4', 'C2xC2']
    for l in [p for p in range(3, lmax) if all(p % q for q in range(2, int(p**0.5)+1)) and p % 4 == 3]:
        h = int(pari.qfbclassno(-l))
        for n in range(1, nmax):
            base = h * (l - 1)//2 * l**delta(l, n)
            if l % 8 == 3 and d % base == 0: G.append('C%d' % l**n)
            if l % 8 == 3 and not (l == 3 and n == 1) and d % (3*base) == 0: G.append('C%d' % (2*l**n))
            if l == 3 and n == 1: G.append('C6')
            if l % 8 == 7 and d % base == 0: G.append('C%d' % (2*l**n))
    return sorted(set(G), key=lambda s: (('x' in s), int(s.split('C')[-1])))
for d in [1, 3, 5, 7, 9]:
    print(d, groups(d))
