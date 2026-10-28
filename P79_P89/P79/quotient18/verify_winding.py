#!/usr/bin/env python3
"""Read-only exact rational replay of the quotient's finite-image criterion."""
from pathlib import Path
import re,json,sys
import sympy as S
H=Path(__file__).resolve().parents[2]/'data/P79/quotient18';s=(H/'finite_image_winding_complete_data.m').read_text()
def expr(n):
 t=re.search(re.escape(n)+r':=(.*?);',s,re.S);assert t,n
 return S.sympify(t[1].replace('true','True').replace('false','False'))
def mat(n):
 t=re.search(re.escape(n)+r':=Matrix\(Rationals\(\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert t,n
 r,c=map(int,t.group(1,2));a=S.sympify(t[3]);assert len(a)==r*c;return S.Matrix(r,c,a)
A,D,AM=map(mat,('Hecke3','Diamond3','MSHecke3'));e=S.Matrix(1,36,expr('MSWindingProjection'));ff=expr('T3Factors');flags=list(expr('WindingZeroFlags'));mins=expr('DiamondMinimalPolynomials');x=S.symbols('x')
polys=[S.Poly(sum(a*x**i for i,a in enumerate(f)),x)for f in ff];p=S.Poly(S.prod(f.as_expr()for f in polys),x);assert A.charpoly(x).as_expr()==p.as_expr()
assert D**3==S.eye(18)and D!=S.eye(18)and(D-S.eye(18)).rank()==12
replayed=[]
for f,mp,flag in zip(polys,mins,flags):
 evalf=S.zeros(18)
 for a in f.all_coeffs():evalf=evalf*A+a*S.eye(18)
 null=evalf.T.nullspace();B=S.Matrix.vstack(*[v.T for v in null]);assert B.rows==f.degree()
 md=S.zeros(18)
 for a in reversed(mp):md=md*D+a*S.eye(18)
 assert B*md==S.zeros(B.rows,18)
 q=S.exquo(p,f);v=S.zeros(1,36)
 for a in q.all_coeffs():v=v*AM+a*e
 zero=v==S.zeros(1,36);replayed.append(bool(zero));assert bool(zero)==bool(flag)
 if zero:assert B*(D-S.eye(18))==S.zeros(B.rows,18)
assert replayed==[True,False,False]
result={'status':'PASS','genus':18,'T3_factor_degrees':[f.degree()for f in polys],'winding_zero_flags':replayed,'diamond3_minimal_polynomials':[[str(a)for a in v]for v in mins],'diamond3_invariant_dimension':6,'alpha_rank':12,'alpha_annihilates_every_winding_zero_factor':True,'finite_image_operator':'diamond3-1'}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
print('NONIC79_FINITE_IMAGE_WINDING_INDEPENDENT_PASS',json.dumps(result,sort_keys=True))
