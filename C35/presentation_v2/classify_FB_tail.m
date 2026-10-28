/* Every degree<=9 place is already in the PROVEN degree9 factor basis. */
AttachSpec("mdmagma/v2/mdmagma.spec");
load "common/pair_af_place_linear.m";
x,y:=Explode(coords);rr:=(x^2*y-x*y+y-1)/(x*(x*y-1));ss:=(x*y-y+1)/(x*y);bb:=rr*ss*(rr-1);cc:=ss*(rr-1);
EE:=EllipticCurve([1-cc,-bb,-bb,F!0,F!0]);PP:=EE![0,0];assert 35*PP eq EE!0 and 5*PP ne EE!0 and 7*PP ne EE!0;
jF:=jInvariant(EE);discF:=Discriminant(EE);cs:=Poles(jF);
assert {*Degree(P):P in cs*} eq {*1^^12,4^^3,6^^2,12*};assert &+[-Valuation(jF,P)*Degree(P):P in cs] eq 576;
function Key35(P)
 key:=[Degree(P),IsFinite(P) select 1 else 0];
 if IsFinite(P) then return key cat [Integers()!a:a in Coefficients(Minimum(P))];end if;return key;
end function;
buckets:=AssociativeArray();
for i in [1..m] do key:=Key35(fb[i]);if IsDefined(buckets,key) then Append(~buckets[key],i);else buckets[key]:=[i];end if;end for;
function FBIndex35(P)
 key:=Key35(P);if not IsDefined(buckets,key) then return 0;end if;
 for i in buckets[key] do if fb[i] eq P then return i;end if;end for;return 0;
end function;
assert &and[FBIndex35(fb[i]) eq i:i in [1..m]];
c0:=Support(Divisor([P:P in cs|Degree(P) eq 1][1]))[1];ic0:=FBIndex35(c0);assert ic0 gt 0;
ci:=[FBIndex35(P):P in cs];smallcuspindices:=[i:i in ci|i gt 0];assert #smallcuspindices eq 17;
qc:=[q(Z.i-Degree(fb[i])*Z.ic0):i in smallcuspindices];HQ:=sub<G|qc>;
printf "NONIC35_FB_CUSP_IMAGE partial_invariants=%o missing_cusp_degrees=%o time=%o\n",Invariants(HQ),[Degree(cs[i]):i in [1..#cs]|ci[i] eq 0],Cputime(tt);
if Invariants(HQ) ne [13,109148520] then
 missing:=[P:P in cs|FBIndex35(P) eq 0];assert #missing eq 1 and Degree(missing[1]) eq 12;
 P12:=missing[1];e12:=Valuation(discF,P12);assert Abs(e12) eq 1;
 ds,es:=Support(Divisor(discF));assert Seqset(ds) subset Seqset(cs);
 cl12:=G!0;
 for i in [1..#ds] do if ds[i] eq P12 then continue;end if;j:=FBIndex35(ds[i]);assert j gt 0;cl12-:=e12*es[i]*q(Z.j-Degree(ds[i])*Z.ic0);end for;
 Append(~qc,cl12);HQ:=sub<G|qc>;
end if;
assert Invariants(HQ) eq [13,109148520] and 10240*#HQ eq KnownOrder;
printf "NONIC35_FB_FULL_RATIONAL_IMAGE invariants=%o index=10240 proven=true time=%o\n",Invariants(HQ),Cputime(tt);
cuspSet:=Seqset(cs);ncindices:=[[i:i in [1..m]|Degree(fb[i]) eq d and fb[i] notin cuspSet]:d in [1..9]];
assert [#a:a in ncindices] eq [0,0,12,0,36,126,348,813,2240];
classes:=[[q(Z.i-d*Z.ic0):i in ncindices[d]]:d in [1..9]];
assert &and[&and[Eltseq(a)[6] eq 0:a in row]:row in classes];
printf "NONIC35_ALL_SMALL_PLACE_CLASSES_FROM_FB counts=%o no_class_map_search=true time=%o\n",[#a:a in classes],Cputime(tt);
load "C35/known_octic35.m";
NN,fc,xc,yc:=Explode(KnownOctic35);RQ<zq>:=PolynomialRing(Rationals());f8:=Polynomial(GF(3),RQ!fc);assert Degree(f8) eq 8 and IsIrreducible(f8);K8<a8>:=ext<GF(3)|f8>;
xy8:=[Evaluate(Polynomial(GF(3),RQ!co),a8):co in [xc,yc]];xx,yy:=Explode(xy8);r8:=(xx^2*yy-xx*yy+yy-1)/(xx*(xx*yy-1));s8:=(xx*yy-yy+1)/(xx*yy);b8:=r8*s8*(r8-1);c8:=s8*(r8-1);
E8:=EllipticCurve([1-c8,-b8,-b8,K8!0,K8!0]);P8:=E8![0,0];assert 35*P8 eq E8!0 and 5*P8 ne E8!0 and 7*P8 ne E8!0;
knownkeys:={};knownclasses:={G|};knownFB:=[];
function PairKey35(v)return Minimum([Eltseq(v[1]^(3^i),GF(3)) cat Eltseq(v[2]^(3^i),GF(3)):i in [0..7]]);end function;
for k in [1..17] do
 if GCD(k,35) ne 1 then continue;end if;v:=TateNormalForm_xy(E8,k*P8);key:=PairKey35(v);if key in knownkeys then continue;end if;Include(~knownkeys,key);
 Q:=C48PairAFPlaceLinear(coords,v);assert Degree(Q) eq 8;ii:=FBIndex35(Q);assert ii gt 0 and ii in ncindices[8];Append(~knownFB,ii);Include(~knownclasses,q(Z.ii-8*Z.ic0));
end for;
assert #knownkeys eq 6 and #knownclasses eq 6 and #Seqset(knownFB) eq 6;
printf "NONIC35_KNOWN_OCTIC_FB_ORBIT places=6 classes=6 time=%o\n",Cputime(tt);
out:=Open("C35/data/presentation_v2/class_data_fb.json","w");
fprintf out,"{\"relations\":%o,\"rational_cusp_generators\":%o,\"known_octic_classes\":%o,\"known_octic_FB_indices\":%o,\"frobenius_coefficients\":%o,\"noncuspidal_counts\":%o,\"tests\":[",[Eltseq(r):r in Rows(RelationMatrix(G))],[Eltseq(a):a in qc],[Eltseq(a):a in knownclasses],knownFB,FrobeniusData,[#a:a in classes];
first:=true;tested:=0;unexpected:=[];good8:={G|};
procedure Test35(typ,inds,a,~first,~tested,~unexpected,~good8)
 d:=&+typ;inH:=a in HQ;isKnown:=d eq 8 and a in knownclasses;
 if not first then fprintf out,",";end if;first:=false;
 fprintf out,"{\"degrees\":%o,\"indices\":%o,\"class\":%o,\"in_H\":%o,\"known_octic\":%o}",typ,inds,Eltseq(a),inH select "true" else "false",isKnown select "true" else "false";tested+:=1;
 if inH then if d eq 8 then Include(~good8,a);if not isKnown then Append(~unexpected,<typ,inds>);end if;else Append(~unexpected,<typ,inds>);end if;end if;
end procedure;
for d in [3,5,6,7,8,9] do
 for i in [1..#classes[d]] do Test35([d],[i],classes[d][i],~first,~tested,~unexpected,~good8);end for;
 printf "NONIC35_FB_PURE_TESTS degree=%o tested=%o unexpected=%o time=%o\n",d,#classes[d],#unexpected,Cputime(tt);
end for;
for i in [1..#classes[3]],j in [i..#classes[3]] do Test35([3,3],[i,j],classes[3][i]+classes[3][j],~first,~tested,~unexpected,~good8);end for;
for i in [1..#classes[5]],j in [1..#classes[3]] do Test35([5,3],[i,j],classes[5][i]+classes[3][j],~first,~tested,~unexpected,~good8);end for;
for i in [1..#classes[6]],j in [1..#classes[3]] do Test35([6,3],[i,j],classes[6][i]+classes[3][j],~first,~tested,~unexpected,~good8);end for;
for i in [1..#classes[3]],j in [i..#classes[3]],k in [j..#classes[3]] do Test35([3,3,3],[i,j,k],classes[3][i]+classes[3][j]+classes[3][k],~first,~tested,~unexpected,~good8);end for;
fprintf out,"],\"complete\":true}\n";delete out;
assert tested eq 5961;
printf "NONIC35_FB_D9_CLASSIFICATION_DONE tests=%o unexpected=%o known_octic_classes=%o complete=%o time=%o\n",tested,#unexpected,#good8,#unexpected eq 0 and good8 eq knownclasses,Cputime(tt);
quit;
