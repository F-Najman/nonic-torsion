N:=73;field_poly:=[ 1, 1, 0, 0, 0, 0, 0, 1 ];SourceJPolynomials:=[* *];TargetJPolynomials:=[* *];SourceECoefficients:=[* *];
Append(~SourceJPolynomials,[ 1, 0, 0, 1, 1, 1, 0, 1 ]);Append(~TargetJPolynomials,[ 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1 ]);Append(~SourceECoefficients,[
[ 1, 0, 0, 0, 0, 0, 0 ],
[ 1, 0, 0, 0, 1, 0, 0 ],
[ 0, 0, 0, 0, 0, 0, 0 ],
[ 0, 0, 0, 0, 0, 0, 0 ],
[ 0, 1, 0, 0, 0, 0, 1 ]
]);
print "NONIC73_NATIVE3_TARGET_DATA_READY";
