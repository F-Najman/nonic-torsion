#!/usr/bin/env python3
"""Read-only exact rational replay of the quotient's finite-image criterion."""
from pathlib import Path
import re,json,sys
import sympy as S
H=Path(__file__).resolve().parents[1]/'data'/'quotient19';s=(H/'finite_image_winding_data.m').read_text()
def expr(n):
 t=re.search(re.escape(n)+r':=(.*?);',s,re.S);assert t,n
 return S.sympify(t[1].replace('true','True').replace('false','False'))
def mat(n):
 t=re.search(re.escape(n)+r':=Matrix\(Rationals\(\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert t,n
 r,c=map(int,t.group(1,2));a=S.sympify(t[3]);assert len(a)==r*c;return S.Matrix(r,c,a)
A,D,AM=map(mat,('HeckeSplit','Diamond5','MSHeckeSplit'));e=S.Matrix(1,38,expr('MSWindingProjection'));ff=expr('HeckeFactors');flags=list(expr('WindingZeroFlags'));mins=expr('DiamondMinimalPolynomials');x=S.symbols('x')
polys=[S.Poly(sum(a*x**i for i,a in enumerate(f)),x)for f in ff];p=S.Poly(S.prod(f.as_expr()for f in polys),x);assert A.charpoly(x).as_expr()==p.as_expr()
assert D**4==S.eye(19)and D**2!=S.eye(19)and(D-S.eye(19)).rank()==14
replayed=[]
for f,mp,flag in zip(polys,mins,flags):
 evalf=S.zeros(19)
 for a in f.all_coeffs():evalf=evalf*A+a*S.eye(19)
 null=evalf.T.nullspace();B=S.Matrix.vstack(*[v.T for v in null]);assert B.rows==f.degree()
 md=S.zeros(19)
 for a in reversed(mp):md=md*D+a*S.eye(19)
 assert B*md==S.zeros(B.rows,19)
 q=S.exquo(p,f);v=S.zeros(1,38)
 for a in q.all_coeffs():v=v*AM+a*e
 zero=v==S.zeros(1,38);replayed.append(bool(zero));assert bool(zero)==bool(flag)
 if zero:assert B*(D-S.eye(19))==S.zeros(B.rows,19)
assert replayed==[False,False,False,True,False,False] and expr('AlphaExponent')==1
result={'status':'PASS','genus':19,'separating_Hecke_factor_degrees':[f.degree()for f in polys],'winding_zero_flags':replayed,'diamond5_minimal_polynomials':[[str(a)for a in v]for v in mins],'diamond5_invariant_dimension':5,'alpha_rank':14,'alpha_annihilates_every_winding_zero_factor':True,'finite_image_operator':'diamond5-1'}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
print('NONIC73_FINITE_IMAGE_WINDING_INDEPENDENT_PASS',json.dumps(result,sort_keys=True))
