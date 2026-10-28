#!/usr/bin/env python3
"""Independent exact Sage/FLINT replay of the finite-image winding criterion."""
from pathlib import Path
import ast,re,json,sys,time
from sage.all import QQ,matrix,PolynomialRing
H=Path(__file__).resolve().parents[1]/'data'/'quotient45';t=time.monotonic();s=(H/'finite_image_winding_data.m').read_text()
def numerical(n):
 if isinstance(n,ast.Constant):assert type(n.value)in(int,bool);return n.value
 if isinstance(n,(ast.List,ast.Tuple)):return[numerical(a)for a in n.elts]
 if isinstance(n,ast.UnaryOp):
  assert isinstance(n.op,(ast.UAdd,ast.USub));return numerical(n.operand)*(-1 if isinstance(n.op,ast.USub)else 1)
 if isinstance(n,ast.BinOp):assert isinstance(n.op,ast.Div);return QQ(numerical(n.left))/QQ(numerical(n.right))
 raise AssertionError(type(n))
def expr(name):
 v=re.search(re.escape(name)+r':=(.*?);',s,re.S);assert v,name
 return numerical(ast.parse(v[1].replace('true','True').replace('false','False'),mode='eval').body)
def mat(name):
 v=re.search(re.escape(name)+r':=Matrix\(Rationals\(\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert v,name
 return matrix(QQ,int(v[1]),int(v[2]),numerical(ast.parse(v[3],mode='eval').body))
A,D,AM=[mat(n)for n in['HeckeSplit','Diamond2','MSHeckeSplit']];assert A.nrows()==45 and AM.nrows()==90
e=matrix(QQ,1,90,expr('MSWindingProjection'));R=PolynomialRing(QQ,'x');ff=[R(a)for a in expr('HeckeFactors')];p=R(1)
for f in ff:assert f.is_irreducible();p*=f
assert A.charpoly()==p and AM.charpoly()==p**2 and p.is_squarefree();I=A.parent().one();assert D**11==I and D!=I and(D-I).rank()==40 and A*D==D*A
flags=expr('WindingZeroFlags');replayed=[];positive=[]
for f,flag in zip(ff,flags):
 q=p//f;v=matrix(QQ,1,90)
 for a in reversed(q.list()):v=v*AM+a*e
 zero=v.is_zero();assert zero==flag;replayed.append(bool(zero))
 if zero:
  B=A.parent().zero()
  for a in reversed(q.list()):B=B*A+a*I
  assert B.rank()==f.degree()and((D-I)*B).is_zero();positive.append(int(f.degree()))
assert positive==[2]and expr('AlphaExponent')==1
out={'status':'PASS','genus':45,'factor_degrees':[int(f.degree())for f in ff],'winding_zero_flags':replayed,'winding_zero_factors_fixed_by_diamond':positive,'alpha_rank':40,'finite_image_operator':'diamond2-1','elapsed_seconds':time.monotonic()-t}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print('NONIC67H3_FINITE_IMAGE_WINDING_REPLAY_PASS',json.dumps(out,sort_keys=True))
