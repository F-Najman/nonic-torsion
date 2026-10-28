SetColumns(0);SetQuitOnError(true);SetMemoryLimit(8*10^9);tt:=Cputime();
load "C44/data/canonical_data.m";
P:=ProjectiveSpace(GF(3),35);R:=CoordinateRing(P);mons:=[R.(p[1])*R.(p[2]):p in pairs];eqs:=[&+[v[i]*mons[i]:i in [1..666]]:v in Rows(QuadRelations)];IC:=ideal<R|eqs>;GB:=GroebnerBasis(IC);LM:=[LeadingMonomial(f):f in GB];
function StdMons(dd)
 ss:={R!1};for k in [1..dd] do cand:={m*R.i:m in ss,i in [1..36]};ss:={m:m in cand|not exists{l:l in LM|IsDivisibleBy(m,l)}};assert #ss eq (k eq 1 select 36 else (2*k-1)*35);end for;return Sort(Setseq(ss));
end function;
S2:=StdMons(2);S4:=StdMons(4);assert #S2 eq 105 and #S4 eq 245;SS4:={@m:m in S4@};
out:=Open("C44/data/multiplication_22_data.m","w");fprintf out,"Std2Exponents:=%o;Std4Exponents:=%o;Tensor22Blocks:=[* *];\n",[Exponents(m):m in S2],[Exponents(m):m in S4];Flush(out);saved:=AssociativeArray();
for i in [1..105] do
 rows:=[];
 for j in [1..105] do
  if j lt i then rr:=saved[<j,i>];else
   nf:=NormalForm(S2[i]*S2[j],IC);rr:=[GF(3)!0:k in [1..245]];mm:=Monomials(nf);cc:=Coefficients(nf);
   for k in [1..#mm] do ii:=Index(SS4,mm[k]);assert ii gt 0;rr[ii]:=cc[k];end for;saved[<i,j>]:=rr;
  end if;Append(~rows,rr);
 end for;
 MM:=Matrix(GF(3),rows);fprintf out,"Append(~Tensor22Blocks,Matrix(GF(3),105,245,%o));\n",[Integers()!x:x in Eltseq(MM)];Flush(out);
 if i mod 10 eq 0 then printf "NONIC44_RING22_BLOCK index=%o time=%o\n",i,Cputime(tt);end if;
end for;
fprintf out,"assert #Tensor22Blocks eq 105;\n";delete out;printf "NONIC44_RING22_COMPLETE R2=105 R4=245 products11025=true symmetric_products5565=true time=%o\n",Cputime(tt);quit;
