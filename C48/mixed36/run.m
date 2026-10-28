SetColumns(0);SetQuitOnError(true);SetMemoryLimit(6*10^9);tt:=Cputime();
AttachSpec("mdmagma/v2/mdmagma.spec");load "C48/degree8/two_isogeny_48_to_2_24.m";load "C48/data/mixed36/inputs.m";
Z:=FreeAbelianGroup(4);Q,qq:=quo<Z|[Z!v:v in QRelations]>;assert #Q eq 4608;actions:=[hom<Q->Q|[qq(Z!v):v in images]>:images in ActionImages];assert &and[#Kernel(g) eq 1 and IsSurjective(g):g in actions];
R:=PolynomialRing(GF(5));K<z>:=ext<GF(5)|R!Source6Polynomial>;
function Key(bc)return Minimum([Eltseq(bc[1]^(5^i),GF(5)) cat Eltseq(bc[2]^(5^i),GF(5)):i in [0..5]]);end function;
classes:=AssociativeArray();sources:=[* *];repclasses:=[];keys:={};
for i in [1..342] do
 bc:=[Evaluate(R!v,z):v in Source6Representatives[i]];b,c:=Explode(bc);E:=EllipticCurve([1-c,-b,-b,K!0,K!0]);P:=E![0,0];assert 48*P eq E!0 and 24*P ne E!0 and 16*P ne E!0;ET,L:=C48TwoIsogenousModuli(E,P);bcT:=DerickxNormalForm_bc(ET,L`P,L`Q);
 canon:=Canonical6Keys[i][2..13];candidates:=[];
 for j in [1..16] do
  eps,eta,a:=Explode(ActionLabels[j]);cur:=DerickxNormalForm_bc(ET,L`P+eps*(12*L`Q),a*L`Q+eta*L`P);
  if Key(cur) eq canon then Append(~candidates,(qq(Z!Canonical6Classes[i]))@@actions[j]);end if;
 end for;assert #candidates gt 0 and #Seqset(candidates) eq 1;cl:=candidates[1];Append(~repclasses,cl);
 for a in [1,5,7,11,13,17,19,23] do
  sbc:=TateNormalForm_bc(E,a*P);ky:=Key(sbc);ta:=Min(a mod 24,24-(a mod 24));j:=Index(ActionLabels,[0,0,ta]);assert j gt 0;v:=actions[j](cl);
  if IsDefined(classes,ky) then assert classes[ky] eq v;else classes[ky]:=v;Include(~keys,ky);Append(~sources,[*i,a,[[Integers()!u:u in Eltseq(t,GF(5))]:t in sbc],Eltseq(v)*]);end if;
 end for;
end for;assert #sources eq 2712;
passes:=[];tested:=0;for i in [1..2],j in [1..#sources] do
 v:=qq(Z!Source3RepClasses[i])+qq(Z!sources[j][4]);tested+:=1;if v eq Q!0 then Append(~passes,[i,j]);end if;
end for;assert tested eq 5424;
out:=Open("C48/data/mixed36/results.py","w");fprintf out,"field_polynomial=%o\nsource6_representative_classes=%o\nsource6_all_rows=%o\npassing_cases=%o\ntested=5424\ncomplete=True\n",Source6Polynomial,[Eltseq(v):v in repclasses],sources,passes;delete out;
printf "NONIC48_MIXED36_DONE allsource6=2712 cases5424=true survivors=%o allsourceclasses_consistent=true complete=true time=%o\n",#passes,Cputime(tt);quit;
