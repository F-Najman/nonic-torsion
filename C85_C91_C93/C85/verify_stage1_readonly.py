#!/usr/bin/env python3
"""Independent exact Tate/unit/finite-recovery replay; no native launch."""
from pathlib import Path
import ast,re,json,time,sys,hashlib
from sage.all import GF,PolynomialRing,LaurentSeriesRing,EllipticCurve,matrix,vector
H=Path(__file__).resolve().parents[1]/'data/C85';t0=time.monotonic();s=(H/'canonical_unit_stage1_data.m').read_text()
def lit(n):return ast.literal_eval(re.search(re.escape(n)+r':=(.*?);',s,re.S)[1])
def mat(n):
 m=re.search(re.escape(n)+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,n;r,c=map(int,m.group(1,2));a=ast.literal_eval(m[3]);assert len(a)==r*c;return[a[i*c:(i+1)*c]for i in range(r)]
def bits(a):return sum(int(v)<<i for i,v in enumerate(a))
def spanrank(rows):
 p={}
 for a in rows:
  while a:
   k=a.bit_length()-1
   if k in p:a^=p[k]
   else:p[k]=a;break
 return len(p)
def combine(v,rows):
 a=0
 while v:
  b=v&-v;a^=rows[b.bit_length()-1];v^=b
 return a
def mul(a,b,mask):
 if a.bit_count()<b.bit_count():a,b=b,a
 z=0
 while b:
  k=b&-b;z^=a<<(k.bit_length()-1);b^=k
 return z&mask
qrows=mat('qcoeffs');Q=[bits(a)for a in qrows];assert len(Q)==225 and all(len(a)==1300 for a in qrows);assert spanrank([a&((1<<577)-1)for a in Q])==225
D=matrix(GF(2),mat('Diamond2Matrix'));assert D**8==D.parent().one()and D**4!=D.parent().one()and(D-D.parent().one()).rank()==196
D3=matrix(GF(2),mat('Diamond3Matrix'));assert D3**16==D3.parent().one()and D3**8!=D3.parent().one()and D3*D==D*D3
LS=LaurentSeriesRing(GF(2),'q',default_prec=1800);q=LS.gen();xp=q/(1+q)**2;yp=q**2/(1+q)**3;s1=LS(0);a6=LS(0)
for n in range(1,1800//85+2):
 if 85*n+1<1800:
  z=q**(85*n+1);xp+=z/(1+z)**2;yp+=z**2/(1+z)**3
 if 85*n-1<1800:
  w=q**(85*n-1);xp+=w/(1+w)**2;yp+=w/(1+w)**3
 if 85*n<1800:
  z=q**(85*n);s1+=GF(2)(n)*z/(1+z);a6+=GF(2)((5*n**3+7*n**5)//12)*z/(1+z)
yp+=s1;assert(yp**2+xp*yp+xp**3+s1*xp+a6).valuation()>=1800

def units(ainv,x,y):
 a1,a2,a3,a4,a6=ainv;aa3=2*y+a3+a1*x;aa2=3*x+a2;aa4=3*x*x+2*x*a2+a4-a1*y;aaa1=2*aa4/aa3+a1;aaa2=aa2-(aa4/aa3)**2-a1*aa4/aa3;b=-aaa2**3/aa3**2;c=-(aaa1*aaa2-aa3)/aa3;r=b/c;ss=c*c/(b-c);t=r*ss-2*r+1;return((ss-r)/t,t/(ss*ss-ss-r+1))
xu,yu=units([LS(1),LS(0),LS(0),s1,a6],xp,yp);assert xu.valuation()==yu.valuation()==-1 and xu.prec()>=1300 and yu.prec()>=1300
low=lit('UnitSeriesLow');assert low==0
assert lit('UnitXSeries')==[int(xu[n])for n in range(1300)]and lit('UnitYSeries')==[int(yu[n])for n in range(1300)]
RX,RY=mat('UnitXRelations'),mat('UnitYRelations');assert len(RX)==136 and len(RY)==131 and all(len(v)==450 for v in RX+RY)
mask=(1<<1299)-1
for u,rels in[(xu,RX),(yu,RY)]:
 uu=bits([int(u[n])for n in range(-1,1300)]);products=[mul(uu,a,(1<<1300)-1)>>1 for a in Q]
 for v in rels:assert(combine(bits(v[:225]),Q)^combine(bits(v[225:]),products))&mask==0
 assert spanrank([bits(v)for v in rels])==len(rels)
# Bidegrees are an upper bound for pole degrees on the normalized component.
model=H.parents[2]/'mdmagma/models_X1_n/FFFc85.txt';text=''.join(model.read_text().split());terms={}
for term in re.findall(r'[+-]?[^+-]+',text):
 term=term.lstrip('+-').rstrip(';');c=1;i=j=0
 for a in term.split('*'):
  if a.startswith('x'):i=int(a.split('^')[1])if'^'in a else 1
  elif a.startswith('y'):j=int(a.split('^')[1])if'^'in a else 1
  else:c*=int(a)
 terms[i,j]=terms.get((i,j),0)^(c&1)
terms={k:v for k,v in terms.items()if v};assert(max(i for i,j in terms),max(j for i,j in terms))==(95,90);assert 1298>2*225-2+95
R=PolynomialRing(GF(2),'a');mp=R(lit('SampleFieldPolynomial'));assert mp.is_irreducible()and mp.degree()==9;F=GF(512,'a',modulus=mp)
def elem(v):return F(R(v))
ainv=list(map(elem,lit('SampleEllipticCoefficients')));E=EllipticCurve(F,ainv);P=E(list(map(elem,lit('SampleLevelPoint'))));assert P!=E(0)and 85*P==E(0)and 5*P!=E(0)and 17*P!=E(0)and E.cardinality()==510
xy=units(ainv,P[0],P[1]);assert xy==tuple(map(elem,lit('SampleUnitCoordinates')))
def recover(xy):
 x,y=xy;A=matrix(F,[[F(v[i])+x*v[225+i]for i in range(225)]for v in RX]+[[F(v[i])+y*v[225+i]for i in range(225)]for v in RY]);assert A.rank()==224;z=A.right_kernel().basis()[0];a=next(t for t in z if t);return vector(F,[t/a for t in z])
z=recover(xy);assert z==vector(F,list(map(elem,lit('SampleCanonicalCoordinates'))));P2=2*P;z2=recover(units(ainv,P2[0],P2[1]));dz=D.change_ring(F)*z;dz/=next(a for a in dz if a);diz=(D**7).change_ring(F)*z;diz/=next(a for a in diz if a);signs=[e for e,v in[(1,dz),(-1,diz)]if v==z2];assert len(signs)==1
# Certify the diamond convention globally, beyond a single finite sample.
# The Tate point u=q^2 is exactly2P on E(q^85). Unit relations at D^{-1}Q
# have zero order beyond deg(K)+max unit pole degree, hence hold globally.
k=2;xp2=q**k/(1+q**k)**2;yp2=q**(2*k)/(1+q**k)**3+s1
for n in range(1,1800//85+2):
 if 85*n+k<1800:
  zq=q**(85*n+k);xp2+=zq/(1+zq)**2;yp2+=zq**2/(1+zq)**3
 if 85*n-k<1800:
  wq=q**(85*n-k);xp2+=wq/(1+wq)**2;yp2+=wq/(1+wq)**3
assert(yp2**2+xp2*yp2+xp2**3+s1*xp2+a6).valuation()>=1800
u2=units([LS(1),LS(0),LS(0),s1,a6],xp2,yp2);DQ=[combine(bits(row),Q)for row in(D**7).rows()];diamond_orders=[]
for u,rels in zip(u2,[RX,RY]):
 v=int(u.valuation());assert v<=0;pr=min(1300+v,int(u.prec())+1);assert pr-1>543
 uc=bits([int(u[n])for n in range(v,pr)]);mm=(1<<(pr-v))-1;pp=[mul(uc,a,mm)for a in DQ]
 for row in rels:assert((combine(bits(row[:225]),DQ)<<(-v))^combine(bits(row[225:]),pp))&mm==0
 diamond_orders.append(pr-1)
k=3;xp3=q**k/(1+q**k)**2;yp3=q**(2*k)/(1+q**k)**3+s1
for n in range(1,1800//85+2):
 if 85*n+k<1800:
  zq=q**(85*n+k);xp3+=zq/(1+zq)**2;yp3+=zq**2/(1+zq)**3
 if 85*n-k<1800:
  wq=q**(85*n-k);xp3+=wq/(1+wq)**2;yp3+=wq/(1+wq)**3
u3=units([LS(1),LS(0),LS(0),s1,a6],xp3,yp3);DQ3=[combine(bits(row),Q)for row in(D3**15).rows()];diamond3_orders=[]
for u,rels in zip(u3,[RX,RY]):
 v=int(u.valuation());assert v<=0;pr=min(1300+v,int(u.prec())+1);assert pr-1>543
 uc=bits([int(u[n])for n in range(v,pr)]);mm=(1<<(pr-v))-1;pp=[mul(uc,a,mm)for a in DQ3]
 for row in rels:assert((combine(bits(row[:225]),DQ3)<<(-v))^combine(bits(row[225:]),pp))&mm==0
 diamond3_orders.append(pr-1)
out={'status':'PASS','canonical_genus':225,'q_precision':1300,'Sturm_weight2':576,'canonical_basis_rank':225,'unit_relations':[136,131],'global_relation_section_degree_bounds':[538,543],'global_zero_order_lower_bound':1298,'sample_finite_rank':224,'sample_level_order':85,'diamond_order':8,'diamond_minus_identity_rank':196,'canonical_matrix_corresponds_to_point_diamond2_exponent':signs[0],'global_diamond_point2_certificate_zero_orders':diamond_orders,'global_diamond_point3_certificate_zero_orders':diamond3_orders,'diamond3_order':16,'native_data_sha256':hashlib.sha256(s.encode()).hexdigest(),'elapsed_seconds':time.monotonic()-t0}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print('NONIC85_CANONICAL_UNIT_STAGE1_INDEPENDENT_PASS',json.dumps(out,sort_keys=True))
