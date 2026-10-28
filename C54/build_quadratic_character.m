SetColumns(0);SetQuitOnError(true);tt:=Cputime();
load "C54/preflight_ordinary_module9.m";
FF5:=FunctionField(C5);AF5,am5:=AlgorithmicFunctionField(FF5);AF5:=RationalExtensionRepresentation(AF5);coords:=[AF5!am5(FF5.i):i in [1..2]];
G5AF,_,phi:=ClassGroup(AF5:Proof:=true);assert G5AF cmpeq G5;
load "C54/quadratic_character.m";
printf "NONIC54_QUADRATIC_BUILD_DONE complete=true time=%o\n",Cputime(tt);quit;
