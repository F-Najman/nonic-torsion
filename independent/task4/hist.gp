\\ Weighted histogram of group structures E(F_q) ~ Z/n1 x Z/n2 over all F_q-isomorphism classes of
\\ elliptic curves, weight 1/#Aut_{F_q}(E).  q = p^e.  Independent enumeration (own code).
\\ Inputs (set before read): p, e, K (number of chunks), k (chunk index 0..K-1), outfile, watch (vector of N), detfile
\\ (parisizemax set on the command line)
q = p^e; g = ffgen(q, 'a);
pw = vector(e, i, g^(i-1));
H = Map();
addw(n1, n2, w) = { my(key = [n1, n2], old); if(mapisdefined(H, key, &old), mapput(H, key, old + w), mapput(H, key, w)); };
structof(E) = { my(G = ellgroup(E)); if(#G == 0, [1,1], if(#G == 1, [1, G[1]], [G[2], G[1]])); };
det(E, j, w, tag) = { my(s = structof(E)); addw(s[1], s[2], w);
  if(#watch, for(i = 1, #watch, if(s[2] % watch[i] == 0,
     write(detfile, [p, e, tag, watch[i], s, w, if(type(j)=="t_FFELT", poldegree(minpoly(j)), 1), Str(j)])))); };
elt(idx) = { my(d = digits(idx, p), n = #d); sum(i = 1, n, d[n + 1 - i] * pw[i]) + 0*g; };
{
\\ nonsquare (odd p) and trace-one element (p = 2)
if(p > 2, ns = 0; until(!issquare(ns) && ns != 0, ns = random(g)));
if(p == 2, delta = 0; until(trace(delta) == 1, delta = random(g)));
lo = (q * k) \ K; hi = (q * (k + 1)) \ K;
mass = 0;
for(idx = lo, hi - 1,
  j = elt(idx);
  if(j == 0 || j == 1728*g^0, next);
  if(p == 2,
    E1 = ellinit([1, 0, 0, 0, 1/j]); E2 = ellinit([1, delta, 0, 0, 1/j]);
  ,
    E0 = ellinit(ellfromj(j));
    b2 = E0.b2; b4 = E0.b4; b6 = E0.b6;
    E1 = ellinit([0, b2/4, 0, b4/2, b6/4]);
    E2 = ellinit([0, ns*b2/4, 0, ns^2*b4/2, ns^3*b6/4]);
  );
  if(E1.j != j || E2.j != j, error("j mismatch"));
  if(ellcard(E1) + ellcard(E2) != 2*q + 2, error("twist mismatch"));
  det(E1, j, 1/2, "gen"); det(E2, j, 1/2, "gen"); mass += 1;
);
\\ special j-invariants: handled in chunk 0 only
if(k == 0,
  if(p == 2,
    \\ j = 0: y^2 + a3 y = x^3 + a4 x + a6, slice a3 in F_q^*/cubes; weight 1/(#R3 q^2)
    r = gcd(3, q - 1); pr = ffprimroot(g); R = vector(r, i, pr^(i-1));
    w = 1/(r * q^2); m0 = 0;
    for(i = 1, r, for(i4 = 0, q - 1, a4 = elt(i4); for(i6 = 0, q - 1, a6 = elt(i6);
       E = ellinit([0, 0, R[i], a4, a6]); if(E.j != 0, error("j0")); det(E, 0, w, "j0"); m0 += w)));
    mass += m0; print("p=2 j=0 slice mass ", m0);
  );
  if(p == 3,
    r = gcd(4, q - 1); pr = ffprimroot(g); R = vector(r, i, pr^(i-1));
    w = 1/(r * q); m0 = 0;
    for(i = 1, r, for(i6 = 0, q - 1, a6 = elt(i6);
       E = ellinit([0, 0, 0, R[i], a6]); if(E.j != 0, error("j0")); det(E, 0, w, "j0"); m0 += w));
    mass += m0; print("p=3 j=0 slice mass ", m0);
  );
  if(p >= 5,
    pr = ffprimroot(g);
    r = gcd(6, q - 1); m0 = 0;
    for(i = 1, r, E = ellinit([0, 0, 0, 0, pr^(i-1)]); if(E.j != 0, error("j0")); det(E, 0, 1/r, "j0"); m0 += 1/r);
    r = gcd(4, q - 1);
    for(i = 1, r, E = ellinit([0, 0, 0, pr^(i-1), 0]); if(E.j != 1728, error("j1728")); det(E, 1728, 1/r, "j1728"); m0 += 1/r);
    mass += m0; print("p>=5 special mass ", m0, " (expect 2)");
  );
);
L = Mat(H);
write(outfile, ["MASS", p, e, k, K, mass]);
for(i = 1, matsize(L)[1], write(outfile, [p, e, k, L[i,1][1], L[i,1][2], L[i,2]]));
write(outfile, ["CHUNK_DONE", p, e, k, K]);
print("done p=", p, " e=", e, " chunk ", k, "/", K, " mass ", mass);
}
quit;
