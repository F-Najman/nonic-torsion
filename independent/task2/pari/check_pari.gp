\\ PARI check: exact torsion via elltors over the nonic field; Codex (p,f,h) triples via ellap at primes of residue degree f.
read("vh_pari.gp");
{
for(i = 1, #ents,
  my(en = ents[i], k = en[1], N = en[3], ex = en[4], exy = en[5], tri = en[6]);
  my(f = subst(ex, y, 0)); f = f / pollead(f, x);
  my(d = 1); for(j = 0, poldegree(f), d = lcm(d, denominator(polcoef(f, j))));
  my(g = d^9 * subst(f, x, x/d));   \\ monic integral, root = d*a
  g = subst(g, x, 'z);
  my(nf = nfinit(g));
  my(a = Mod('z / d, g));
  my(cy = polcoef(exy, 1, y), c0 = polcoef(exy, 0, y));
  my(x0 = a, y0 = -subst(c0, x, a) / subst(cy, x, a));
  my(r = (x0^2*y0 - x0*y0 + y0 - 1)/(x0^2*y0 - x0), s = (x0*y0 - y0 + 1)/(x0*y0));
  my(b = r*s*(r-1), c = s*(r-1));
  my(E = ellinit([1-c, -b, -b, 0, 0], nf));
  my(T = elltors(E));
  my(P = [0, 0]); my(ok = ellisoncurve(E, P) && ellorder(E, P) == N);
  my(res = []);
  for(t = 1, #tri, my(p = tri[t][1], f0 = tri[t][2], h = tri[t][3], dec = idealprimedec(nf, p), hits = []);
     for(m = 1, #dec, my(pr = dec[m]); if(pr.f == f0 && pr.e == 1,
        my(hh = if(ellglobalred(E)[1] % p == 0 && idealval(nf, ellglobalred(E)[1], pr) > 0, "bad", pr.p^pr.f + 1 - ellap(E, pr)));
        hits = concat(hits, [hh])));
     res = concat(res, [[p, f0, h, hits]]));
  print("PARI entry ", k, " N=", N, " elltors=", T[2], " P_order_ok=", ok, " codex_triples: ", res);
);
}
\\ explicit C2 x C26 curve (Prop 3.44)
{
my(g = z^9+5*z^8+15*z^7+14*z^6-14*z^5-25*z^4-7*z^3+5*z^2+4*z+1, nf = nfinit(g), a = Mod(z, g));
my(b = (62*a^8+1386*a^7+6639*a^6+18950*a^5+20738*a^4-6786*a^3-25017*a^2-14822*a-2411)/907);
my(c = (-3864*a^8-22538*a^7-75946*a^6-112918*a^5-25630*a^4+94470*a^3+104296*a^2+46530*a+8797)/907);
my(E = ellinit([0, c, 0, b*(1-b-c), 0], nf));
my(T = elltors(E));
print("PARI C2xC26 elltors=", T[2], " ord(0,0)=", ellorder(E,[0,0]), " ord(b,b)=", ellorder(E,[b,b]));
my(res = []);
foreach(idealprimedec(nf, 7), pr, res = concat(res, [[pr.f, pr.e, pr.p^pr.f + 1 - ellap(E, pr)]]));
print("  primes above 7 [f,e,#E]: ", res);
res = [];
foreach(idealprimedec(nf, 3), pr, res = concat(res, [[pr.f, pr.e, pr.p^pr.f + 1 - ellap(E, pr)]]));
print("  primes above 3 [f,e,#E]: ", res);
print("  poly mod 3 factorization: ", factormod(g, 3)[,1]~);
print("  poly mod 7 factorization: ", factormod(g, 7)[,1]~);
}
print("PARI_DONE");
quit;
