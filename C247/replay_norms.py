#!/usr/bin/env python3
"""Independent exact determinant/norm parity replay over Q."""
from pathlib import Path
import json
from fractions import Fraction
from sympy import Matrix,Rational
P=Path(__file__).resolve().parent
expected={(6,3):(44,235561454207215818721),(3,3):(44,63860140169462891961),(2,17):(20,8744990250917207615976649)}
for (o,q),(d,n) in expected.items():
 D=json.loads((P/'data'/f'norm_order{o}_q{q}.json').read_text())
 assert D['order']==o and D['q']==q and D['rational_matrix_dimension']==d
 assert len(D['entries'])==d*d
 E=Matrix(d,d,[Rational(a,b) for a,b in D['entries']]);det=E.det(method='domain-ge')
 assert det==n==D['exact_determinant_norm'] and n%2==1
 print(f'NONIC247_INDEPENDENT_NORM_PASS order={o} q={q} rational_dimension={d} determinant={n} odd=true')
print('NONIC247_NORM_PARITY_REPLAY_COMPLETE')
