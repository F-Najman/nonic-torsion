from pathlib import Path
from fractions import Fraction
import re,json,sys
from sage.all import QQ,GF,PolynomialRing,NumberField,EllipticCurve,prime_range,version
P=Path(__file__).resolve().parents[1]/'data'/'explicit'
s=(P/'example.txt').read_text().replace('[*','[').replace('*]',']');s=re.sub(r'(-?\d+)/(\d+)',r'Fraction(\1,\2)',s);data=eval(s,{'__builtins__':{},'Fraction':Fraction});pol,uk,bc=data
R=PolynomialRing(QQ,'z');zz=R.gen();f=R(pol);assert f.degree()==9
R3=PolynomialRing(GF(3),'z');assert R3(f).is_irreducible()
K=NumberField(f,'a');a=K.gen();b,c=[sum(QQ(x)*a**i for i,x in enumerate(row)) for row in bc];E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P2=E(0,0);Q26=E(b,b)
assert 2*P2==E(0) and P2!=E(0) and 26*Q26==E(0) and 13*Q26!=E(0) and 2*Q26!=E(0) and P2!=13*Q26
print('C226_SAGE_EXACT_LEVEL_PASS degree=9 full2=true primitive26=true',flush=True)
reductions=[];obstructions={};factor7=[];order3=None
for p in prime_range(3,102):
 Rp=PolynomialRing(GF(p),'z');ff=Rp(f)
 if ff.gcd(ff.derivative()).degree():continue
 for h,m in ff.factor():
  assert m==1;d=h.degree()
  if d==1:L=GF(p);z=-h[0]/h[1]
  else:L=GF(p**d,'w',modulus=h);z=L.gen()
  bb,cc=[sum(L(x.numerator)*L(x.denominator)**(-1)*z**i for i,x in enumerate(row)) for row in bc]
  EE=EllipticCurve(L,[0,cc,0,bb*(1-bb-cc),0])
  if not EE.discriminant():continue
  pp=EE(0,0);q=EE(bb,bb);roots=[L(0),(13*q)[0],(pp+13*q)[0]]
  for i,r in enumerate(roots):
   other=[roots[j] for j in range(3) if j!=i];A=cc+3*r;B=(r-other[0])*(r-other[1]);delta=other[0]-other[1]
   assert delta**2==A*A-4*B
   for j,x in enumerate([delta,-delta]):
    val=2*x*(delta**2-A*x)
    if val and val**((L.cardinality()-1)//2)==-L(1) and (i,j) not in obstructions:
     obstructions[i,j]=dict(prime=int(p),modulus=[int(t) for t in h.list()],value=[int(t) for t in val.polynomial().list()] if d>1 else [int(val)])
  if p==3:
   order3=int(EE.cardinality());assert order3==19604
  if p==7:
   n=int(EE.cardinality());assert (d==2 and n==52) or (d==3 and n==364);factor7.append(d);reductions.append(dict(degree=int(d),order=n))
 if p>=7 and len(obstructions)==6:break
assert sorted(factor7)==[2,2,2,3] and len(obstructions)==6
out=dict(degree=9,full_torsion=[2,26],prime3_order=order3,reduction7=reductions,four_division_obstructions={str(i+1)+','+str(j+1):v for (i,j),v in sorted(obstructions.items())},complete=True)
if '--readonly' in sys.argv:assert json.loads((P/'independent_certificate.json').read_text())==out
else:(P/'independent_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C226_SAGE_FULL_TORSION_AND_ALL_ISOGENY_LOCAL_OBSTRUCTIONS_PASS',json.dumps(out),flush=True)
