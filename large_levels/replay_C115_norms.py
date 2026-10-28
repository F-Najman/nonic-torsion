#!/usr/bin/env python3
"""Independent exact resultant/norm/parity verification for the C115 proof."""
import json
from pathlib import Path
import sympy as s
H=Path(__file__).resolve().parent/'data';t=s.symbols('t');out=[]
for order in [44,11,22]:
 d=json.loads((H/f'C115_norm_order{order}.json').read_text())
 assert d['N']==115 and d['q']==3 and d['relative_dimension']==10
 f=sum(s.Rational(c)*t**i for i,c in enumerate(d['field_polynomial']))
 v=sum(s.Rational(c)*t**i for i,c in enumerate(d['determinant_coefficients']))
 assert s.degree(f,t)==d['coefficient_degree']==20 and s.Poly(f,t).is_irreducible
 assert f==s.cyclotomic_poly(44,t)
 exact=s.resultant(f,v,t);assert exact==int(d['norm']) and exact%2==1
 assert d['norm_power']==20//s.totient(order)
 out.append({'character_order':order,'norm':str(exact),'norm_power':d['norm_power'],'odd':True})
assert sum(10*int(s.totient(o)) for o in [44,11,22])==400
(H/'C115_norm_replay.json').write_text(json.dumps({'degree':9,'N':115,'image_dimension':400,'verified':True,'rows':out},indent=2)+'\n')
print('NONIC_C115_NORM_RESULTANT_REPLAY_PASS image_dimension=400 characters=3 all_norms_odd=true complete=true')
