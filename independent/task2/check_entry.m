// Independent verification of a degree-9 torsion realization.
// Usage: magma -b Idx:=k check_entry.m   (k = 1..19: van Hoeij degree-9 entries; k = 20: explicit C2xC26 curve)
SetColumns(0);
SetMemoryLimit(24*10^9);
if not assigned Idx then Idx := 1; elif Type(Idx) eq MonStgElt then Idx := StringToInteger(Idx); end if;
R<x,y> := PolynomialRing(Rationals(), 2);
Px<t> := PolynomialRing(Rationals());
load "independent/task2/vh_deg9.m";
t0 := Cputime();
if Idx le #VH then
  e := VH[Idx]; line := e[1]; N := e[2]; extra := e[3]; eqx := e[4]; eqxy := e[5];
  f := Px ! UnivariatePolynomial(eqx);
  assert Degree(f) eq 9 and IsIrreducible(f);
  K<a> := NumberField(f);
  assert Degree(eqxy, y) eq 1;
  cy := Coefficient(eqxy, y, 1); c0 := Coefficient(eqxy, y, 0);
  ev := func<g | Evaluate(g, [a, K!0])>;
  x0 := a; y0 := -ev(c0)/ev(cy);
  assert Evaluate(eqxy, [x0, y0]) eq 0;
  r := (x0^2*y0 - x0*y0 + y0 - 1)/(x0^2*y0 - x0);
  s := (x0*y0 - y0 + 1)/(x0*y0);
  b := r*s*(r-1); c := s*(r-1);
  E := EllipticCurve([1-c, -b, -b, 0, 0]);
  Pt := E![0,0];
  gens := [Pt]; claimed := [N];
  printf "ENTRY idx=%o srcline=%o N=%o extra=%o\n", Idx, line, N, extra;
else
  N := 52;
  f := t^9+5*t^8+15*t^7+14*t^6-14*t^5-25*t^4-7*t^3+5*t^2+4*t+1;
  assert IsIrreducible(f);
  K<a> := NumberField(f);
  b := (62*a^8+1386*a^7+6639*a^6+18950*a^5+20738*a^4-6786*a^3-25017*a^2-14822*a-2411)/907;
  c := (-3864*a^8-22538*a^7-75946*a^6-112918*a^5-25630*a^4+94470*a^3+104296*a^2+46530*a+8797)/907;
  E := EllipticCurve([0, c, 0, b*(1-b-c), 0]);
  P2 := E![0,0]; Q := E![b,b];
  gens := [P2, Q];
  printf "ENTRY idx=%o C2xC26 explicit (Prop 3.44)\n", Idx;
end if;
printf "FIELDPOLY %o\n", f;
assert Discriminant(E) ne 0;
// exact orders of the given points
for Pt in gens do
  n := Order(Pt);
  printf "POINT_ORDER %o\n", n;
end for;
if Idx le #VH then
  ok := IsZero(N*gens[1]) and &and[ not IsZero((N div q)*gens[1]) : q in PrimeDivisors(N)];
  printf "ORDER_CHECK N=%o exact=%o\n", N, ok;
else
  Q := gens[2]; P2 := gens[1];
  ok := IsZero(26*Q) and not IsZero(13*Q) and not IsZero(2*Q) and IsZero(2*P2) and not IsZero(P2) and P2 ne 13*Q;
  printf "ORDER_CHECK C2xC26 exact=%o\n", ok;
end if;
// --- field data
OK := MaximalOrder(K);
dK := Discriminant(OK);
printf "FIELD_DISC %o = %o\n", dK, Factorization(dK);
r1, r2 := Signature(K); printf "SIGNATURE (%o,%o)\n", r1, r2;
G, roots, data := GaloisGroup(K);
printf "GALOIS order=%o transitive_id=9T%o\n", #G, TransitiveGroupIdentification(G);
subs := Subfields(K);
printf "SUBFIELD_DEGREES %o\n", [Degree(S[1]) : S in subs];
for S in subs do if Degree(S[1]) eq 3 then printf "CUBIC_SUBFIELD %o disc %o\n", DefiningPolynomial(S[1]), Discriminant(MaximalOrder(S[1])); end if; end for;
// --- j-invariant
j := jInvariant(E);
mj := MinimalPolynomial(j);
printf "J_DEGREE %o\n", Degree(mj);
if Degree(mj) le 3 then printf "J_MINPOLY %o\n", mj; end if;
printf "J_INTEGRAL %o\n", IsIntegral(j);
if Degree(mj) eq 1 then printf "J_RATIONAL %o\n", Rationals()!j; end if;
// CM: CM j-invariants are algebraic integers
if not IsIntegral(j) then
  printf "CM false (j not integral)\n";
else
  cm := HasComplexMultiplication(E);
  if cm then _, D := HasComplexMultiplication(E); printf "CM true D=%o\n", D; else printf "CM false (HasComplexMultiplication)\n"; end if;
end if;
// --- base change from a proper subfield (degree 1 or 3)
c4E, c6E := Explode(cInvariants(E));
if Degree(mj) eq 1 then
  jQ := Rationals()!j;
  assert jQ ne 0 and jQ ne 1728;
  E0 := EllipticCurve([0,0,0,-27*jQ*(jQ-1728), -54*jQ*(jQ-1728)^2]);
  assert jInvariant(E0) eq jQ;
  c40, c60 := Explode(cInvariants(E0));
  d := c6E*c40/(c4E*c60);
  nd := Norm(d);
  bc := IsSquare(d/nd);
  printf "BASECHANGE subfield_deg=1 j_in_subfield=true twist_class_in_subfield=%o\n", bc;
  if bc then
    E0t := MinimalModel(QuadraticTwist(E0, nd));
    printf "BASECHANGE_ISOMORPHIC %o\n", IsIsomorphic(ChangeRing(E0t, K), E);
    printf "BASECHANGE_MODEL_OVER_Q %o conductor %o\n", aInvariants(E0t), Conductor(E0t);
    printf "BASECHANGE_TORSION_OVER_Q %o\n", Invariants(TorsionSubgroup(E0t));
  end if;
end if;
for S in subs do
  F := S[1]; emb := S[2];
  if Degree(F) eq 9 or Degree(F) eq 1 then continue; end if;
  // is j in emb(F)?
  // robust test: j in emb(F) iff minimal polynomial of j has a root in F mapping to j
  found := false;
  for rt in Roots(ChangeRing(mj, F)) do
    if emb(rt[1]) eq j then found := true; jF := rt[1]; end if;
  end for;
  if not found then printf "BASECHANGE subfield_deg=%o j_in_subfield=false\n", Degree(F); continue; end if;
  if jF eq 0 or jF eq 1728 then printf "BASECHANGE subfield_deg=%o j=%o special, skipped\n", Degree(F), jF; continue; end if;
  E0 := EllipticCurve([F!0,0,0, -27*jF*(jF-1728)/1, -54*jF*(jF-1728)^2]); // j(E0)=jF
  assert jInvariant(E0) eq jF;
  c40, c60 := Explode(cInvariants(E0));
  d := c6E*emb(c40)/(c4E*emb(c60));   // E ~ twist of E0 by d (mod squares)
  if Degree(F) eq 1 then
    issq := IsSquare(d);
    // over Q: d in Q* K*^2 ?  Q*/Q*^2 -> K*/K*^2 ; test N_{K/Q}(d) * d is a square (odd degree 9)
    nd := Norm(d);
    bc := IsSquare(d/nd);
    printf "BASECHANGE subfield_deg=1 j_in_subfield=true twist_class_in_subfield=%o\n", bc;
  else
    L := RelativeField(F, K);
    dL := L!d;
    nd := Norm(dL);  // in F
    bc := IsSquare(d/emb(nd));
    printf "BASECHANGE subfield_deg=3 j_in_subfield=true twist_class_in_subfield=%o\n", bc;
  end if;
  if bc then
    // exhibit E0' over F with E0' x K isomorphic to E: E0' = twist of E0 by nd
    E0t := QuadraticTwist(E0, nd);
    E0tK := EllipticCurve([emb(ai) : ai in aInvariants(E0t)]);
    printf "BASECHANGE_ISOMORPHIC %o\n", IsIsomorphic(E0tK, E);
    printf "BASECHANGE_MODEL_OVER_SUBFIELD %o\n", aInvariants(E0t);
    T0, m0 := TorsionSubgroup(E0t);
    printf "BASECHANGE_TORSION_OVER_SUBFIELD %o\n", Invariants(T0);
  end if;
end for;
// --- independent reduction bound: odd primes p not dividing disc(f_int), unramified, Dedekind factorization
fi := f; den := LCM([Denominator(cc) : cc in Coefficients(f)]);
// integral monic model: a' = den*a has minimal poly  den^9 f(t/den)
fint := Px ! (den^9 * Evaluate(f, t/den));
assert IsMonic(fint) and &and[IsIntegral(cc) : cc in Coefficients(fint)];
Zt<u> := PolynomialRing(Integers());
fZ := Zt ! fint; disc := Discriminant(fZ);
ainv := aInvariants(E);
// express a-invariants as polynomials in a' = den*a
coeffpols := [ Px ! Eltseq(ai) : ai in ainv ];
coeffpols := [ Evaluate(cp, t/den) : cp in coeffpols ];
g := 0; used := [];
for p in PrimesInInterval(3, 400) do
  if disc mod p eq 0 then continue; end if;
  if &or[ Denominator(cc) mod p eq 0 : cc in Coefficients(cp), cp in coeffpols ] then continue; end if;
  Fp := GF(p); Fpt<w> := PolynomialRing(Fp);
  for fac in Factorization(Fpt ! fZ) do
    k := Degree(fac[1]); if p^k gt 10^12 then continue; end if;
    Fq := ext<Fp | fac[1]>; rt := Roots(Fpt!fac[1], Fq)[1][1];
    ai := [ Evaluate(ChangeRing(Px!cp, Fp), rt) : cp in coeffpols ];
    okE, Eb := IsEllipticCurve(ai);
    if not okE then continue; end if;   // model singular mod this prime: skip
    h := #Eb; g := GCD(g, h); Append(~used, <p, k, h>);
  end for;
  if #used ge 40 then break; end if;
end for;
printf "REDUCTION_GCD %o from %o places; first: %o\n", g, #used, used[1..Min(8,#used)];
// --- Magma TorsionSubgroup
T, mT := TorsionSubgroup(E);
printf "TORSION_SUBGROUP %o\n", Invariants(T);
printf "DONE idx=%o time=%o\n", Idx, Cputime(t0);
quit;
