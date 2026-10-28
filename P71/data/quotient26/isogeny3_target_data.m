N:=71;SourceDegrees:=[];SourceJPolynomials:=[* *];TargetJPolynomials:=[* *];SourceTargetIndices:=[* *];SourceECoefficients:=[* *];SourceFieldPolynomials:=[* *];
Append(~TargetJPolynomials,[ 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1 ]);
Append(~SourceDegrees,7);Append(~SourceJPolynomials,[ 1, 1, 1, 1, 1, 1, 0, 1 ]);Append(~SourceTargetIndices,[ 1 ]);Append(~SourceECoefficients,[
[ 1, 0, 0, 0, 0, 0, 0 ],
[ 1, 1, 1, 1, 0, 0, 0 ],
[ 0, 0, 0, 0, 0, 0, 0 ],
[ 0, 0, 0, 0, 0, 0, 0 ],
[ 1, 1, 0, 0, 0, 0, 0 ]
]);Append(~SourceFieldPolynomials,[ 1, 1, 0, 0, 0, 0, 0, 1 ]);
Append(~TargetJPolynomials,[ 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1 ]);
Append(~TargetJPolynomials,[ 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1 ]);
Append(~SourceDegrees,8);Append(~SourceJPolynomials,[ 1, 1, 0, 1, 1, 1, 1, 0, 1 ]);Append(~SourceTargetIndices,[ 2, 3 ]);Append(~SourceECoefficients,[
[ 1, 0, 0, 0, 0, 0, 0, 0 ],
[ 0, 0, 0, 0, 0, 0, 0, 0 ],
[ 0, 0, 0, 0, 0, 0, 0, 0 ],
[ 0, 0, 0, 0, 0, 0, 0, 0 ],
[ 0, 1, 0, 1, 1, 1, 0, 0 ]
]);Append(~SourceFieldPolynomials,[ 1, 0, 1, 1, 1, 0, 0, 0, 1 ]);
print "NONIC71_NATIVE3_TARGET_DATA_READY";
