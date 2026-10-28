#!/usr/bin/env python3
"""Direct Tate-unit orders at every cusp of an odd squarefree level."""
from pathlib import Path
import json,time,sys,math,resource,os
resource.setrlimit(resource.RLIMIT_CPU,(120,125));resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));os.nice(10)
os.environ["OPENBLAS_NUM_THREADS"]="1";os.environ["OMP_NUM_THREADS"]="1"
import faulthandler;faulthandler.dump_traceback_later(10,repeat=True)
print("SAGE_IMPORT_START",flush=True)
from sage.all import GF,LaurentSeriesRing,divisors,euler_phi
print("SAGE_IMPORT_COMPLETE",flush=True)
N=int(sys.argv[1])if len(sys.argv)>1 else 93;assert N in[65,91,93];H=Path(__file__).resolve().parents[1]/'data/N93_source_units';t0=time.monotonic();PREC=300

def units(x,y,a4):
 w=(x*x+a4+y)/x;aa=x+w*w+w;b=aa**3/x**2;c=(aa+x)/x;r=b/c;s=c*c/(b+c);tt=r*s+1;return((s+r)/tt,tt/(s*s+s+r+1))
def order2(m):
 if m==1:return 1
 n=1;t=2%m
 while t!=1:n+=1;t=2*t%m
 return n
import faulthandler
faulthandler.dump_traceback_later(10,repeat=True)
records=[]
for d in reversed(divisors(N)):
 d=int(d)
 if d==1:continue
 print("WIDTH_START",d,flush=True);m=N//d;df=order2(m);F=GF(2**df,'a')if df>1 else GF(2);LS=LaurentSeriesRing(F,'q',default_prec=PREC);q=LS.gen();a4=LS(0);a6=LS(0)
 for n in range(1,PREC//d+1):
  z=q**(d*n);a4+=F(n)*z/(1+z);a6+=F((5*n**3+7*n**5)//12)*z/(1+z)
 roots=[];seen=set()
 for b in range(1,m+1):
  if math.gcd(b,m)!=1 or b%m in seen:continue
  roots.append(b);a=b%m
  for _ in range(df):seen.add(a);a=2*a%m
 assert len(seen)==int(euler_phi(m))
 zeta=F.multiplicative_generator()**((F.cardinality()-1)//m)if m>1 else F(1)
 print("FAMILY_READY",d,"roots",roots,"df",df,flush=True)
 for a in range(1,d//2+1):
  if math.gcd(a,d)!=1:continue
  for exponent in roots:
   zz=zeta**exponent;uk=zz*q**a;xp=uk/(1+uk)**2;yp=uk**2/(1+uk)**3+a4
   for n in range(1,PREC//d+2):
    if d*n+a<PREC:z=zz*q**(d*n+a);xp+=z/(1+z)**2;yp+=z*z/(1+z)**3
    if d*n-a<PREC:w=zz**(-1)*q**(d*n-a);xp+=w/(1+w)**2;yp+=w/(1+w)**3
   assert(yp*yp+xp*yp+xp**3+a4*xp+a6).valuation()>=PREC
   x,y=units(xp,yp,a4);assert x.prec()>100 and y.prec()>100;records.append({'width':d,'component':a,'root_exponent_representative':exponent,'residue_degree':df,'orders':[int(x.valuation()),int(y.valuation())]})
print("ALL_NON_WIDTH1_DONE",flush=True)
F=GF(2**order2(N),'a');zeta=F.multiplicative_generator()**((F.cardinality()-1)//N);assert zeta.multiplicative_order()==N
for b in range(1,(N+1)//2):
 if math.gcd(b,N)!=1:continue
 z=zeta**b;xp=z/(1+z)**2;yp=z*z/(1+z)**3;x,y=units(xp,yp,F(0));assert x and y
px=sum(r['residue_degree']*max(0,-r['orders'][0])for r in records);py=sum(r['residue_degree']*max(0,-r['orders'][1])for r in records);width1=int(euler_phi(N))//2;geometric=sum(r['residue_degree']for r in records)+width1
assert geometric==sum(int(euler_phi(d)*euler_phi(N//d))for d in divisors(N))//2
out={'status':'PASS','N':N,'precision':PREC,'pole_degrees':[px,py],'geometric_cusps':geometric,'non_width1_records':records,'width1_geometric_cusps':width1,'width1_orders':[0,0],'plane_model_used':False,'canonical_relations_used':False,'elapsed_seconds':time.monotonic()-t0};saved=json.loads((H/'cusp_pole_degrees_verified_v3.json').read_text());assert {k:v for k,v in saved.items()if k!='elapsed_seconds'}=={k:v for k,v in out.items()if k!='elapsed_seconds'};print('DIRECT_CUSP_POLE_AUDIT_PASS',N,[px,py],geometric,'elapsed',time.monotonic()-t0)
