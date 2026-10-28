SetColumns(0);SetQuitOnError(true);SetMemoryLimit(16*10^9);
AttachSpec("mdmagma/v2/mdmagma.spec");tt:=Cputime();
X:=MDX1(37,GF(2));FF:=FunctionField(Curve(X));AF,am:=AlgorithmicFunctionField(FF);assert Genus(AF) eq 40;
u:=am(X`_coordinates_xy[1]);v:=am(X`_coordinates_xy[2]);K:=BaseRing(AF);
printf "NONIC37_AF_COORDINATES degree=%o u_base=%o v_base=%o u_gen=%o v_gen=%o time=%o\n",Degree(AF),u eq AF!K.1,v eq AF!K.1,u eq AF.1,v eq AF.1,Cputime(tt);
E:=EllipticCurve([am(z):z in X`_E]);P:=E![AF!0,AF!0];xy:=TateNormalForm_xy(E,2*P);
if u eq AF!K.1 and v eq AF.1 then sig0:=hom<K->AF|xy[1]>;dg:=hom<AF->AF|sig0,xy[2]>;else error "unexpected AF coordinates";end if;
printf "NONIC37_DIAMOND_MAP_CONSTRUCTED time=%o\n",Cputime(tt);
Om,om:=SpaceOfHolomorphicDifferentials(AF);B:=[om(vv):vv in Basis(Om)];
CM:=Matrix(GF(2),[Eltseq(Cartier(w)@@om):w in B]);eta:=Differential(u);deta:=Differential(dg(u));MD:=Matrix(GF(2),[Eltseq((dg(w/eta)*deta)@@om):w in B]);I:=IdentityMatrix(GF(2),40);
assert MD^18 eq I and MD*CM eq CM*MD;U:=RowSpace(MD-I);assert U eq RowSpace((CM-I)*(MD-I)) and Dimension(U) eq 38;
printf "NONIC37_ALPHA_CARTIER_IMAGE dimension=%o time=%o\n",Dimension(U),Cputime(tt);

load "common/dt_expansions_local.m";
load "P37/data/input_records.m";
function XYKey(xx,yy)
 d:=LCM(Degree(MinimalPolynomial(xx,GF(2))),Degree(MinimalPolynomial(yy,GF(2))));found:=false;
 for bits in [[1,0],[0,1],[1,1]] do zz:=bits[1]*xx+bits[2]*yy;mp:=MinimalPolynomial(zz,GF(2));if Degree(mp) eq d then bb:=bits;found:=true;break;end if;end for;assert found;
 BM:=Matrix(GF(2),[Eltseq(zz^i,GF(2)):i in [0..d-1]]);assert Rank(BM) eq d;
 cs:=[Eltseq(Solution(BM,Vector(GF(2),Eltseq(vv,GF(2))))):vv in [xx,yy]];
 return [d]cat bb cat[Integers()!cc:cc in Eltseq(mp)]cat[Integers()!cc:cc in &cat cs];
end function;
rec:=Inputs[1];assert rec[1] eq 6 and rec[2] eq 1;PF:=PolynomialRing(GF(2));F<zz>:=ext<GF(2)|PF!rec[3]>;
EE:=EllipticCurve([F!PF!vv:vv in rec[4]]);PN:=EE![F!PF!vv:vv in rec[5]];xx,yy:=Explode(TateNormalForm_xy(EE,PN));key:=XYKey(xx,yy);
zs:=Zeros(Evaluate(MinimalPolynomial(xx,GF(2)),u));pp:=[];
for pl in zs do
 if Degree(pl) ne 6 or Valuation(u,pl) lt 0 or Valuation(v,pl) lt 0 then continue;end if;
 if XYKey(Evaluate(u,pl),Evaluate(v,pl)) eq key then Append(~pp,pl);end if;
end for;assert #pp eq 1;P6:=pp[1];
b:=am(X`_coordinates[1]);c:=am(X`_coordinates[2]);j:=(1+c)^12/(b^3*(b+c*(c+1)^3));cs:=Poles(j);assert Degree(j) eq 684 and {*Degree(pl):pl in cs*} eq {*1^^18,18*};rats:=[pl:pl in cs|Degree(pl) eq 1];
places:=[P6]cat rats;ds:=[6]cat[1:i in [1..18]];caps:=[1]cat[3:i in [1..18]];
f:="P37/data/sextic_cusp3_jets.m";
PrintFile(f,Sprintf("N:=37;D:=9;gx:=40;udim:=38;ds:=%o;caps:=%o;basisjets:=AssociativeArray();localjets:=AssociativeArray();",ds,caps):Overwrite:=true);
PrintFile(f,Sprintf("CartierMatrix:=Matrix(GF(2),40,40,%o);DiamondMatrix:=Matrix(GF(2),40,40,%o);ImageBasis:=Matrix(GF(2),38,40,%o);",Eltseq(CM),Eltseq(MD),Eltseq(BasisMatrix(U))));
for i in [1..19] do
 printf "NONIC37_SEXTIC_CUSP_JET_START i=%o degree=%o cap=%o time=%o\n",i,ds[i],caps[i],Cputime(tt);
 AA:=HorizontalJoin(DifferentialExpansionMatrices(B,places[i],caps[i]));assert Nrows(AA) eq 40 and Ncols(AA) eq ds[i]*caps[i];
 PrintFile(f,Sprintf("basisjets[%o]:=Matrix(GF(2),40,%o,%o);",i,Ncols(AA),Eltseq(AA)));
 for m in [1..caps[i]] do
  JJ:=BasisMatrix(U)*Submatrix(AA,1,1,40,m*ds[i]);assert Rank(JJ) eq m*ds[i];
  PrintFile(f,Sprintf("localjets[<%o,%o>]:=Matrix(GF(2),38,%o,%o);",i,m,Ncols(JJ),Eltseq(JJ)));
 end for;
 printf "NONIC37_SEXTIC_CUSP_JET_DONE i=%o time=%o\n",i,Cputime(tt);
end for;
PrintFile(f,"print \"NONIC37_SEXTIC_CUSP3_JET_DATA_READY\";");
VV,vvm:=DifferentialSpace(Divisor(P6)+3*Divisor(rats[1]));VO:=VectorSpace(GF(2),40);W:=sub<VO|[VO!Eltseq(vvm(vv)@@om):vv in Basis(VV)]>;assert 38-Dimension(U meet W) eq 9;
printf "NONIC37_SEXTIC_CUSP3_NATIVE_RR_CONTROL rank=9 complete=true time=%o\n",Cputime(tt);
print "NONIC37_SEXTIC_CUSP3_EXPORT_COMPLETE";quit;
