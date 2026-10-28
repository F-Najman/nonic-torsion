SetColumns(0);SetQuitOnError(true);tt:=Cputime();
load "C133/data/fi3_data_sharp_direct.m";
assert gx eq 35 and D eq 3 and udim eq 8;
DM:=DiamondMatrix;CM:=CartierMatrix;I:=IdentityMatrix(GF(2),35);
assert DM^3 eq I and DM*CM eq CM*DM;T:=DM*CM;
U0:=RowSpace((T-DM)*(I+DM+DM^2)*(T^4+T^2+I));
U3:=RowSpace((T-DM)*(DM+DM^2)*T);U:=U0+U3;
assert Dimension(U0) eq 2 and Dimension(U3) eq 20 and Dimension(U) eq 22;
for i in [1..#ds],m in [1..3 div ds[i]] do
 A:=Submatrix(basisjets[i],1,1,35,m*ds[i]);
 assert Nullspace(A) eq RowSpace(kernelspaces[<i,m>]);
 assert ImageBasis*A eq localjets[<i,m>];
end for;
function WV(ds,d)
 if #ds eq 0 then return d eq 0 select [[]] else [];end if;result:=[];
 for m in [0..d div ds[1]] do for v in WV(ds[2..#ds],d-m*ds[1]) do Append(~result,[m] cat v);end for;end for;return result;
end function;
vs:=WV(ds,3);assert #vs eq 58;bad:=[];
for v in vs do
 A:=ZeroMatrix(GF(2),35,0);
 for i in [1..#ds] do if v[i] gt 0 then A:=HorizontalJoin(A,Submatrix(basisjets[i],1,1,35,v[i]*ds[i]));end if;end for;
 assert Ncols(A) eq 3;rk:=Rank(BasisMatrix(U)*A);if rk lt 3 then Append(~bad,v);end if;
end for;
printf "NONIC133_SPLIT_NATIVE_REPLAY U0=2 U3=20 combined_U=22 exact_RR_blocks=20 all_cusp_configurations=58 exceptions=%o complete=true time=%o\n",#bad,Cputime(tt);
assert #bad eq 0;quit;
