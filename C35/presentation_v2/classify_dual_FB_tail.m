/* Exact primary-image membership in the rational cuspidal subgroup. */
load "C35/data/presentation_v2/dual_characters_data.m";
assert DualBaseIndex35 eq Index(FBDegrees,1);ic0:=DualBaseIndex35;
assert Matrix(Integers(16),RawRelations)*Transpose(DualMatrix16) eq 0 and Matrix(GF(5),RawRelations)*Transpose(DualMatrix5) eq 0;
assert &and[DualMatrix16[i,ic0] eq 0:i in [1..Nrows(DualMatrix16)]] and &and[DualMatrix5[i,ic0] eq 0:i in [1..2]];
AttachSpec("mdmagma/v2/mdmagma.spec");
load "common/pair_af_place_linear.m";
x,y:=Explode(coords);r:=(x^2*y-x*y+y-1)/(x*(x*y-1));s:=(x*y-y+1)/(x*y);b:=r*s*(r-1);c:=s*(r-1);
E:=EllipticCurve([1-c,-b,-b,F!0,F!0]);P0:=E![0,0];assert 35*P0 eq E!0 and 5*P0 ne E!0 and 7*P0 ne E!0;
jF:=jInvariant(E);cs:=Poles(jF);assert {*Degree(P):P in cs*} eq {*1^^12,4^^3,6^^2,12*};assert &+[-Valuation(jF,P)*Degree(P):P in cs] eq 576;
assert fb[ic0] in cs;
function Key35(P)
 key:=[Degree(P),IsFinite(P) select 1 else 0];if IsFinite(P) then return key cat [Integers()!a:a in Coefficients(Minimum(P))];end if;return key;
end function;
buckets:=AssociativeArray();for i in [1..m] do key:=Key35(fb[i]);if IsDefined(buckets,key) then Append(~buckets[key],i);else buckets[key]:=[i];end if;end for;
function FBIndex35(P)
 key:=Key35(P);if not IsDefined(buckets,key) then return 0;end if;for i in buckets[key] do if fb[i] eq P then return i;end if;end for;return 0;
end function;
assert &and[FBIndex35(fb[i]) eq i:i in [1..m]];
ci:=[FBIndex35(P):P in cs];cuspindices:=[i:i in ci|i gt 0];assert #cuspindices eq 17;
r16:=Nrows(DualMatrix16);G16:=AbelianGroup([16:i in [1..r16]]);G5:=AbelianGroup([5,5]);
col16:=[G16![Integers()!DualMatrix16[i,j]:i in [1..r16]]:j in [1..m]];col5:=[G5![Integers()!DualMatrix5[i,j]:i in [1..2]]:j in [1..m]];
I16:=sub<G16|col16>;I5:=sub<G5|col5>;assert #I16 eq 2^14 and #I5 eq 25;
H16:=sub<G16|[col16[i]:i in cuspindices]>;H5:=sub<G5|[col5[i]:i in cuspindices]>;
assert Invariants(H16) eq [8] and Invariants(H5) eq [5];
printf "NONIC35_PRIMARY_RATIONAL_IMAGES full2_order=16384 full5_order=25 cusp2_order=8 cusp5_order=5 quotient_index=10240 exact=true time=%o\n",Cputime(tt);
cuspSet:=Seqset(cs);ncindices:=[[i:i in [1..m]|Degree(fb[i]) eq d and fb[i] notin cuspSet]:d in [1..9]];
assert [#a:a in ncindices] eq [0,0,12,0,36,126,348,813,2240];
load "C35/known_octic35.m";
NN,fc,xc,yc:=Explode(KnownOctic35);RQ<zq>:=PolynomialRing(Rationals());f8:=Polynomial(GF(3),RQ!fc);assert IsIrreducible(f8) and Degree(f8) eq 8;K8<a8>:=ext<GF(3)|f8>;
xy8:=[Evaluate(Polynomial(GF(3),RQ!co),a8):co in [xc,yc]];xx,yy:=Explode(xy8);r8:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));s8:=(xx*yy-yy+1)/(xx*yy);b8:=r8*s8*(r8-1);c8:=s8*(r8-1);
E8:=EllipticCurve([1-c8,-b8,-b8,K8!0,K8!0]);P8:=E8![0,0];assert 35*P8 eq E8!0 and 5*P8 ne E8!0 and 7*P8 ne E8!0;
function PairKey35(v)return Minimum([Eltseq(v[1]^(3^i),GF(3)) cat Eltseq(v[2]^(3^i),GF(3)):i in [0..7]]);end function;
knownkeys:={};knownFB:=[];
for k in [1..17] do
 if GCD(k,35) ne 1 then continue;end if;v:=TateNormalForm_xy(E8,k*P8);key:=PairKey35(v);if key in knownkeys then continue;end if;Include(~knownkeys,key);
 Q:=C48PairAFPlaceLinear(coords,v);assert Degree(Q) eq 8;ii:=FBIndex35(Q);assert ii in ncindices[8];Append(~knownFB,ii);
end for;
assert #knownkeys eq 6 and #Seqset(knownFB) eq 6;
assert &and[col16[i] in H16 and col5[i] in H5:i in knownFB];
printf "NONIC35_PRIMARY_KNOWN_OCTIC places=6 exact_FB_indices=%o time=%o\n",knownFB,Cputime(tt);
out:=Open("C35/data/presentation_v2/dual_classification.json","w");
fprintf out,"{\"level\":35,\"base_FB_index\":%o,\"cusp_FB_indices\":%o,\"noncuspidal_FB_indices\":%o,\"known_octic_FB_indices\":%o,\"full_primary_image_orders\":[16384,25],\"rational_primary_image_orders\":[8,5],\"tests\":[",ic0,cuspindices,ncindices,knownFB;
first:=true;tested:=0;unexpected:=[];good8:={Integers()|};
procedure Test35(typ,inds,a16,a5,~first,~tested,~unexpected,~good8)
 d:=&+typ;inH:=a16 in H16 and a5 in H5;isKnown:=false;
 if typ eq [8] then isKnown:=ncindices[8][inds[1]] in knownFB;end if;
 if not first then fprintf out,",";end if;first:=false;
 fprintf out,"{\"degrees\":%o,\"indices\":%o,\"image16\":%o,\"image5\":%o,\"in_H\":%o,\"known_octic_place\":%o}",typ,inds,Eltseq(a16),Eltseq(a5),inH select "true" else "false",isKnown select "true" else "false";tested+:=1;
 if inH then if isKnown then Include(~good8,ncindices[8][inds[1]]);else Append(~unexpected,<typ,inds>);end if;end if;
end procedure;
for d in [3,5,6,7,8,9] do
 for i in [1..#ncindices[d]] do j:=ncindices[d][i];Test35([d],[i],col16[j],col5[j],~first,~tested,~unexpected,~good8);end for;
 printf "NONIC35_PRIMARY_PURE_DONE degree=%o places=%o unexpected=%o time=%o\n",d,#ncindices[d],#unexpected,Cputime(tt);
end for;
for i in [1..#ncindices[3]],j in [i..#ncindices[3]] do a:=ncindices[3][i];bb:=ncindices[3][j];Test35([3,3],[i,j],col16[a]+col16[bb],col5[a]+col5[bb],~first,~tested,~unexpected,~good8);end for;
for de in [<5,3>,<6,3>] do d,e:=Explode(de);for i in [1..#ncindices[d]],j in [1..#ncindices[e]] do a:=ncindices[d][i];bb:=ncindices[e][j];Test35([d,e],[i,j],col16[a]+col16[bb],col5[a]+col5[bb],~first,~tested,~unexpected,~good8);end for;end for;
for i in [1..#ncindices[3]],j in [i..#ncindices[3]],k in [j..#ncindices[3]] do a:=ncindices[3][i];bb:=ncindices[3][j];cc:=ncindices[3][k];Test35([3,3,3],[i,j,k],col16[a]+col16[bb]+col16[cc],col5[a]+col5[bb]+col5[cc],~first,~tested,~unexpected,~good8);end for;
fprintf out,"],\"complete\":true}\n";delete out;assert tested eq 5961;
printf "NONIC35_PRIMARY_CLASSIFICATION_DONE tests=5961 unexpected=%o known_octic_places=%o complete=%o time=%o\n",#unexpected,#good8,#unexpected eq 0 and good8 eq Seqset(knownFB),Cputime(tt);
quit;
