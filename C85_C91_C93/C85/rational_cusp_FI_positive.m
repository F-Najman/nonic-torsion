// Formal immersion directly in the saturated canonical q-expansion lattice.
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(4*10^9);tt:=Cputime();
load "C85_C91_C93/data/C85/canonical_unit_stage1_data.m";
Q:=qcoeffs;D2:=Diamond2Matrix;D3:=Diamond3Matrix;II:=IdentityMatrix(GF(2),225);prec:=1300;
function CanonicalHecke(q,D)
 len:=(prec-1)div q+1;pre:=Submatrix(Q,1,1,225,len);assert Rank(pre)eq 225;DQ:=D*Q;rr:=Matrix(GF(2),[[Q[i,q*n+1]+(n mod q eq 0 select (GF(2)!q)*DQ[i,n div q+1]else GF(2)!0):n in[0..len-1]]:i in[1..225]]);T:=Solution(pre,rr);assert T*pre eq rr;return T;
end function;
T2:=CanonicalHecke(2,D2);T3:=CanonicalHecke(3,D3);
S:=CuspidalSubspace(ModularSymbolsH(85,[1],2,1));assert Dimension(S)eq 225;
positive:=[A:A in NewformDecomposition(S)|IsZero(WindingElementProjection(A))];assert #positive eq 2 and &and[Dimension(A)eq 2 and Level(AssociatedNewSpace(A))eq 85:A in positive];
labels:=[Min(u,85-u)where u:=2^a*3^b mod 85:b in[0..15],a in[0..1]];dd:=[D2^a*D3^b:b in[0..15],a in[0..1]];assert #Seqset(labels)eq 32 and labels[1]eq 1;
combined:=ZeroMatrix(GF(2),133,0);used:=[];killcofs:=[* *];flagdata:=[* *];complete:=false;
for n in[3,2]do
 polys:=[CharacteristicPolynomial(HeckeOperator(A,n)):A in positive];kp:=LCM(polys);assert Degree(kp)le 4;AA:=MatrixAlgebra(GF(2),225);mt:=Evaluate(ChangeRing(kp,GF(2)),AA!(n eq 3 select T3 else T2))*(T3-II-D3);rr:=[];
 for k in[1..32]do
  fq:=mt*dd[k]*Q;for j in[1..(k eq 1 select 9 else 4)]do Append(~rr,[fq[i,j+1]:i in[1..225]]);end for;
 end for;
 mat:=Matrix(GF(2),rr);assert Nrows(mat)eq 133 and Ncols(mat)eq 225;combined:=HorizontalJoin(combined,mat);Append(~used,n);Append(~killcofs,[Integers()!a:a in Eltseq(kp)]);Append(~flagdata,[[Integers()!c:c in Eltseq(f)]:f in polys]);
 printf "NONIC85_COTANGENT_FI_PROGRESS n=%o killpolynomial=%o operatorrank=%o cumulative_rank=%o rows133=true columns=%o time=%o\n",n,kp,Rank(mt),Rank(combined),Ncols(combined),Cputime(tt);
 if Rank(combined)eq 133 then complete:=true;break;end if;
end for;
root:="C85_C91_C93/data/C85/";
PrintFile(root cat "FI_positive_rows.txt",Sprintf("85 32 133 %o",Ncols(combined)):Overwrite:=true);for row in Rows(combined)do PrintFile(root cat "FI_positive_rows.txt",&cat[IntegerToString(Integers()!a):a in Eltseq(row)]);end for;
PrintFile(root cat "FI_positive_data.m",Sprintf("N:=85;SourceGenus:=225;DiamondLabels:=%o;HeckeIndices:=%o;KillPolynomialCoefficients:=%o;WindingFactorData:=%o;CanonicalHecke2:=Matrix(GF(2),225,225,%o);CanonicalHecke3:=Matrix(GF(2),225,225,%o);Rows:=133;Columns:=%o;Rank:=%o;Complete:=%o;",labels,used,killcofs,flagdata,Eltseq(T2),Eltseq(T3),Ncols(combined),Rank(combined),complete):Overwrite:=true);
if complete then print "NONIC85_COTANGENT_FI_COMPLETE all_normalized_degree9_configs_rank9=true";else print "NONIC85_COTANGENT_FI_PARTIAL no_exclusion_inferred";end if;quit;
