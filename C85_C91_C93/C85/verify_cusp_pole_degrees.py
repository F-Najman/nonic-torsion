#!/usr/bin/env python3
"""Direct exact Tate cusp pole degrees at squarefree level85."""
from pathlib import Path
import json,time,sys,math
from sage.all import GF,LaurentSeriesRing
H=Path(__file__).resolve().parents[1]/'data/C85';t0=time.monotonic();PREC=600

def units(x,y,a4):
 w=(x*x+a4+y)/x;aa=x+w*w+w;b=aa**3/x**2;c=(aa+x)/x;r=b/c;s=c*c/(b+c);tt=r*s+1;return((s+r)/tt,tt/(s*s+s+r+1))
def family(d,F):
 LS=LaurentSeriesRing(F,'q',default_prec=PREC);q=LS.gen();s1=LS(0);a6=LS(0)
 for n in range(1,PREC//d+1):
  z=q**(d*n);s1+=F(n)*z/(1+z);a6+=F((5*n**3+7*n**5)//12)*z/(1+z)
 return LS,q,s1,a6
records=[]
for d,df,roots in[(85,1,[1]),(17,4,[1]),(5,8,[1,3])]:
 F=GF(2**df,'a')if df>1 else GF(2);LS,q,a4,a6=family(d,F);order=85//d;zeta=F.multiplicative_generator()**((F.cardinality()-1)//order)if order>1 else F(1)
 for a in range(1,d//2+1):
  if math.gcd(a,d)!=1:continue
  for exponent in roots:
   zz=zeta**exponent;uk=zz*q**a;xp=uk/(1+uk)**2;yp=uk**2/(1+uk)**3+a4
   for n in range(1,PREC//d+2):
    if d*n+a<PREC:z=zz*q**(d*n+a);xp+=z/(1+z)**2;yp+=z*z/(1+z)**3
    if d*n-a<PREC:w=zz**(-1)*q**(d*n-a);xp+=w/(1+w)**2;yp+=w/(1+w)**3
   assert(yp*yp+xp*yp+xp**3+a4*xp+a6).valuation()>=PREC
   x,y=units(xp,yp,a4);assert x.prec()>100 and y.prec()>100;records.append({'width':d,'component':a,'root_exponent_representative':exponent,'residue_degree':df,'orders':[int(x.valuation()),int(y.valuation())]})
F=GF(256,'a');zeta=F.multiplicative_generator()**3;assert zeta.multiplicative_order()==85
for b in range(1,43):
 if math.gcd(b,85)!=1:continue
 z=zeta**b;xp=z/(1+z)**2;yp=z*z/(1+z)**3;x,y=units(xp,yp,F(0));assert x and y
px=sum(r['residue_degree']*max(0,-r['orders'][0])for r in records);py=sum(r['residue_degree']*max(0,-r['orders'][1])for r in records);geometric=sum(r['residue_degree']for r in records)+32;assert geometric==128 and(px,py)==(90,95)
out={'status':'PASS','N':85,'precision':PREC,'pole_degrees':[px,py],'geometric_cusps':128,'non_width1_records':records,'width1_geometric_cusps':32,'width1_orders':[0,0],'plane_model_used':False,'canonical_relations_used':False,'elapsed_seconds':time.monotonic()-t0};p=H/'cusp_pole_degrees_verified.json'
if '--readonly'in sys.argv:
 old=json.loads(p.read_text());assert{k:v for k,v in old.items()if k!='elapsed_seconds'}=={k:v for k,v in out.items()if k!='elapsed_seconds'}
else:p.write_text(json.dumps(out,indent=2)+'\n')
print('NONIC85_DIRECT_CUSP_POLE_DEGREES_PASS poles90_95 cusps128 elapsed',time.monotonic()-t0,flush=True)
