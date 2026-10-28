/* Small primary character targets; no large integral quotient transformation. */
SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
load "C35/data/presentation_v2/verification_inputs.m";
assert N eq 35 and p eq 3 and CandidateInvariants eq [4,4,8,520,218297040,0];
m:=#FBDegrees;ic0:=Index(FBDegrees,1);assert m eq 3592 and ic0 gt 0;
delete FBBases;delete FieldPolynomialData;delete ModularCoordinatesData;
R16:=Integers(16);A16:=Matrix(R16,RawRelations);unit:=ZeroMatrix(R16,1,m);unit[1,ic0]:=1;A16:=VerticalJoin(A16,unit);
K16:=Nullspace(Transpose(A16));B16:=BasisMatrix(K16);assert #K16 eq 2^14 and A16*Transpose(B16) eq 0;
r16:=Nrows(B16);P16:=AbelianGroup([16:i in [1..r16]]);I16:=sub<P16|[P16![Integers()!B16[i,j]:i in [1..r16]]:j in [1..m]]>;assert #I16 eq 2^14;
printf "NONIC35_DUAL16_READY generators=%o kernel_size=%o image_order=%o exact_relations=true time=%o\n",r16,#K16,#I16,Cputime(tt);
A5:=Matrix(GF(5),RawRelations);unit5:=ZeroMatrix(GF(5),1,m);unit5[1,ic0]:=1;A5:=VerticalJoin(A5,unit5);K5:=Nullspace(Transpose(A5));B5:=BasisMatrix(K5);
assert Dimension(K5) eq 2 and A5*Transpose(B5) eq 0;
P5:=AbelianGroup([5,5]);I5:=sub<P5|[P5![Integers()!B5[i,j]:i in [1..2]]:j in [1..m]]>;assert #I5 eq 25;
out:=Open("C35/data/presentation_v2/dual_characters.json","w");
fprintf out,"{\"level\":35,\"factor_basis_size\":3592,\"base_FB_index\":%o,\"modulus16\":16,\"character_generators16\":%o,\"character_matrix16\":%o,\"image_order16\":16384,\"modulus5\":5,\"character_generators5\":2,\"character_matrix5\":%o,\"image_order5\":25,\"all_relations_annihilated\":true,\"complete\":true}\n",ic0,r16,[[Integers()!B16[i,j]:j in [1..m]]:i in [1..r16]],[[Integers()!B5[i,j]:j in [1..m]]:i in [1..2]];delete out;
out:=Open("C35/data/presentation_v2/dual_characters_data.m","w");fprintf out,"DualBaseIndex35:=%o;\nDualMatrix16:=Matrix(Integers(16),%o,3592,%o);\nDualMatrix5:=Matrix(GF(5),2,3592,%o);\n",ic0,r16,[Integers()!a:a in Eltseq(B16)],[Integers()!a:a in Eltseq(B5)];delete out;
printf "NONIC35_PRIMARY_CHARACTERS_DONE full_primary_image_order=409600 native_proven_Pic_order=14529850982400 complete=true time=%o\n",Cputime(tt);quit;
