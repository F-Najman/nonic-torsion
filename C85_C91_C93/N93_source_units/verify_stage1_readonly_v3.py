#!/usr/bin/env python3
"""Independent exact Tate/unit/finite-recovery replay; no native launch."""
from pathlib import Path
import ast,re,json,time,sys,hashlib,os
os.environ["OPENBLAS_NUM_THREADS"]="1";os.environ["OMP_NUM_THREADS"]="1"
from sage.all import GF,PolynomialRing,LaurentSeriesRing,EllipticCurve,matrix,vector
H=Path(__file__).resolve().parents[1]/'data/N93_source_units';t0=time.monotonic();s=(H/'canonical_unit_data.m').read_text()
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
qrows=mat('qcoeffs');Q=[bits(a)for a in qrows];assert len(Q)==261 and all(len(a)==1800 for a in qrows);assert spanrank([a&((1<<642)-1)for a in Q])==261
D=matrix(GF(2),mat('Diamond1Matrix'));E5=matrix(GF(2),mat('Diamond2Matrix'));assert D**10==E5**6==D.parent().one()and E5**3==D**5
LS=LaurentSeriesRing(GF(2),'q',default_prec=2300);q=LS.gen();xp=q/(1+q)**2;yp=q**2/(1+q)**3;s1=LS(0);a6=LS(0)
for n in range(1,2300//93+2):
 if 93*n+1<2300:
  z=q**(93*n+1);xp+=z/(1+z)**2;yp+=z**2/(1+z)**3
 if 93*n-1<2300:
  w=q**(93*n-1);xp+=w/(1+w)**2;yp+=w/(1+w)**3
 if 93*n<2300:
  z=q**(93*n);s1+=GF(2)(n)*z/(1+z);a6+=GF(2)((5*n**3+7*n**5)//12)*z/(1+z)
yp+=s1;assert(yp**2+xp*yp+xp**3+s1*xp+a6).valuation()>=2300

def units(ainv,x,y):
 a1,a2,a3,a4,a6=ainv;aa3=2*y+a3+a1*x;aa2=3*x+a2;aa4=3*x*x+2*x*a2+a4-a1*y;aaa1=2*aa4/aa3+a1;aaa2=aa2-(aa4/aa3)**2-a1*aa4/aa3;b=-aaa2**3/aa3**2;c=-(aaa1*aaa2-aa3)/aa3;r=b/c;ss=c*c/(b-c);t=r*ss-2*r+1;return((ss-r)/t,t/(ss*ss-ss-r+1))
xu,yu=units([LS(1),LS(0),LS(0),s1,a6],xp,yp);assert xu.valuation()==yu.valuation()==-1 and xu.prec()>=1800 and yu.prec()>=1800
low=lit('UnitSeriesLow');assert low==0
assert lit('UnitXSeries')==[int(xu[n])for n in range(1800)]and lit('UnitYSeries')==[int(yu[n])for n in range(1800)]
RX,RY=mat('UnitXRelations'),mat('UnitYRelations');assert len(RX)==162 and len(RY)==157 and all(len(v)==522 for v in RX+RY)
mask=(1<<1799)-1
for u,rels in[(xu,RX),(yu,RY)]:
 uu=bits([int(u[n])for n in range(-1,1800)]);products=[mul(uu,a,(1<<1800)-1)>>1 for a in Q]
 for v in rels:assert(combine(bits(v[:261]),Q)^combine(bits(v[261:]),products))&mask==0
 assert spanrank([bits(v)for v in rels])==len(rels)
pole=json.loads((H/'cusp_pole_degrees_verified_v3.json').read_text());assert pole['status']=='PASS'and pole['pole_degrees']==[100,105]and pole['geometric_cusps']==120
assert 1798>2*261-2+105
# Certify the diamond convention globally, beyond a single finite sample.
# The Tate point u=q^2 is exactly2P on E(q^93). Unit relations at D^{-1}Q
# have zero order beyond deg(K)+max unit pole degree, hence hold globally.
diamond_orders=[]
for k,DM in[(2,D**9),(5,E5**5)]:
 xp2=q**k/(1+q**k)**2;yp2=q**(2*k)/(1+q**k)**3+s1
 for n in range(1,2300//93+2):
  if 93*n+k<2300:
   zq=q**(93*n+k);xp2+=zq/(1+zq)**2;yp2+=zq**2/(1+zq)**3
  if 93*n-k<2300:
   wq=q**(93*n-k);xp2+=wq/(1+wq)**2;yp2+=wq/(1+wq)**3
 assert(yp2**2+xp2*yp2+xp2**3+s1*xp2+a6).valuation()>=2300
 u2=units([LS(1),LS(0),LS(0),s1,a6],xp2,yp2);DQ=[combine(bits(row),Q)for row in DM.rows()]
 for u,rels in zip(u2,[RX,RY]):
  v=int(u.valuation());assert v<=0;pr=min(1800+v,int(u.prec())+1);assert pr-1>625
  uc=bits([int(u[n])for n in range(v,pr)]);mm=(1<<(pr-v))-1;pp=[mul(uc,a,mm)for a in DQ]
  for row in rels:assert((combine(bits(row[:261]),DQ)<<(-v))^combine(bits(row[261:]),pp))&mm==0
  diamond_orders.append(pr-1)
out={'status':'PASS','canonical_genus':261,'q_precision':1800,'Sturm_weight2':640,'canonical_basis_rank':261,'unit_relations':[162,157],'direct_cusp_pole_bounds':[100,105],'global_relation_section_degree_bounds':[620,625],'global_zero_order_lower_bound':1798,'diamond_orders':[10,6],'global_point_diamonds2_5_inverse':True,'global_diamond_zero_orders':diamond_orders,'native_data_sha256':hashlib.sha256(s.encode()).hexdigest()}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print('NONIC93_CANONICAL_UNIT_STAGE1_INDEPENDENT_PASS',out)
