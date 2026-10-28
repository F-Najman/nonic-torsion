load "common/pair_af_place_linear.m";
function ClassGroupWithAF(C,c0)
 FF:=FunctionField(C);AF,am:=AlgorithmicFunctionField(FF);
 G,_,phi:=ClassGroup(AF:ReductionDivisor:=FunctionFieldDivisor(Divisor(c0)),Proof:=true);
 m2:=map<DivisorGroup(C)->G|D:->phi(FunctionFieldDivisor(D)),g:->CurveDivisor(C,g@@phi)>;
 return G,Inverse(m2),m2,AF,am,phi;
end function;
