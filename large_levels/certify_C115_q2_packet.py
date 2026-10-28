#!/usr/bin/env python3
"""Read-only exact audit of the C115 q2 certificate."""
from pathlib import Path
from fractions import Fraction
import hashlib,json,re,math
import sympy as s
H=Path(__file__).resolve().parent/'data';L=H.parent/'logs';M=json.loads((H/'C115_q2_manifest.json').read_text())
assert M['degree']==9 and M['excluded_cyclic_orders']==[115] and M['complete']
assert M['degree9_principal_bound']==54<Fraction(M['gonality_lower_bound'])
assert pow(16,11,115)==1 and all(pow(16,j,115) not in [1,114] for j in range(1,11))
x=s.symbols('x');dim=0
for o in [44,11,22]:
 d=json.loads((H/f'C115_norm_order{o}.json').read_text());assert d['N']==115 and d['q']==3
 f=sum(s.Rational(c)*x**i for i,c in enumerate(d['field_polynomial']));v=sum(s.Rational(c)*x**i for i,c in enumerate(d['determinant_coefficients']))
 assert f==s.cyclotomic_poly(44,x) and d['norm_power']==20//s.totient(o)
 r=s.resultant(f,v,x);assert r==int(d['norm']) and r%2==1
 dim+=int(s.totient(o))*d['relative_dimension']
assert dim==400
log=(L/'C115_image_even_orders_certified.log').read_text()
assert 'NONIC_115_IMAGE_NO_TWO_PRIMARY_PASS q=3 complete=true' in log
assert log.rstrip().endswith('NONIC_115_IMAGE_EVEN_ORDERS_CERTIFIED_COMPLETED')
assert not re.search('Runtime error|User error|Assertion failed',log)
wind=(L/'inputs/winding115.log').read_text();rr=re.findall(r'WINDING N=115 order=(\d+) conductor=\d+ relative_dimension=(\d+) winding_dimension=(\d+)',wind)
assert sum(int(s.totient(int(o)))*int(d) for o,d,w in rr)==441
assert all(d==w for o,d,w in rr if int(o)>1)
print('NONIC_C115_Q2_READONLY_PASS files='+str(len(M['files']))+' exact_resultants=3 norm_parity=odd image_dimension=400 complete=true')
