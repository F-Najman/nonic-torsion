// Independent analytic-rank check for J1(N) (Kind:="J1") or J1(2,2n) ~ J_H(4n) (Kind:="J22", N:=n).
// For every level M | L (L = N, resp. 4n) and every Galois orbit of even Dirichlet characters chi mod M
// (for J22: conductor(chi) | 2n), decompose the new cuspidal sign+1 modular symbols into newform
// orbits and test whether the winding element has nonzero projection (<=> L(f,1) != 0 for the orbit).
// Weighted dimension sum is compared with the genus of the curve as a completeness check.
SetColumns(0);
SetMemoryLimit(24*10^9);
if not assigned N then N := 27; elif Type(N) eq MonStgElt then N := StringToInteger(N); end if;
if not assigned Kind then Kind := "J1"; end if;
t0 := Cputime();
L := Kind eq "J1" select N else 4*N;
condbound := Kind eq "J1" select N else 2*N;
total := 0; nonfull := [];
for M in Divisors(L) do
  mult := #Divisors(L div M);
  DG := FullDirichletGroup(M);
  reps := GaloisConjugacyRepresentatives(DG);
  for chi0 in reps do
    if IsOdd(chi0) then continue; end if;
    if condbound mod Conductor(chi0) ne 0 then continue; end if;
    chi := MinimalBaseRingCharacter(chi0);
    degchi := Degree(BaseRing(chi)) ;
    MS := ModularSymbols(chi, 2, +1);
    S := NewSubspace(CuspidalSubspace(MS));
    d := Dimension(S);
    if d eq 0 then continue; end if;
    D := NewformDecomposition(S);
    for A in D do
      dA := Dimension(A);
      w := Dimension(WindingSubmodule(A));
      full := (w eq dA);
      total +:= mult*degchi*dA;
      printf "FACTOR M=%o chi_order=%o chi_cond=%o degQchi=%o dim=%o absdim=%o mult=%o winding_full=%o\n",
        M, Order(chi), Conductor(chi), degchi, dA, dA*degchi, mult, full;
      if not full then
        // character values on generators of (Z/LZ)^*, to test diamond-difference claims
        vals := [ <a, chi(a) eq 1> : a in [2..L] | GCD(a, L) eq 1 ];
        Append(~nonfull, <M, Order(chi), Conductor(chi), dA*degchi, [v[1] : v in vals | v[2]]>);
      end if;
    end for;
  end for;
end for;
// genus via the dimension formula over all even characters mod L with conductor | condbound
g := 0;
for chi0 in GaloisConjugacyRepresentatives(FullDirichletGroup(L)) do
  if IsOdd(chi0) or condbound mod Conductor(chi0) ne 0 then continue; end if;
  g +:= Degree(BaseRing(MinimalBaseRingCharacter(chi0))) * DimensionCuspForms(chi0, 2);
end for;
if Kind eq "J1" then printf "GENUS_GAMMA1 %o\n", Genus(Gamma1(N)); end if;
printf "SUMMARY Kind=%o N=%o weighted_dim=%o genus=%o allfull=%o time=%o\n", Kind, N, total, g, #nonfull eq 0, Cputime(t0);
for x in nonfull do
  printf "NONFULL M=%o chi_order=%o chi_cond=%o absdim=%o chi_trivial_at=%o\n", x[1], x[2], x[3], x[4], x[5];
end for;
printf "RANKCHECK_DONE\n";
quit;
