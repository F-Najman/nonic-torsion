SetColumns(0);SetQuitOnError(true);SetMemoryLimit(12*10^9);tt:=Cputime();
load "C54/preflight_ordinary_module9.m";
FF5:=FunctionField(C5);AF5,am5:=AlgorithmicFunctionField(FF5);AF5:=RationalExtensionRepresentation(AF5);coords:=[AF5!am5(FF5.i):i in [1..2]];
G5AF,afinv,phi:=ClassGroup(AF5:Proof:=true);assert G5AF cmpeq G5;A0:=FunctionFieldDivisor(Divisor(c05));assert Degree(A0) eq 1;
