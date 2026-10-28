#!/usr/bin/env python3
"""Independent pole degrees from all Tate cusps; no plane-model dependency."""
from pathlib import Path
import time,json,sys
from sage.all import GF,LaurentSeriesRing,PolynomialRing,QQ,FractionField
H=Path(__file__).resolve().parents[1]/'data/pole_bound_supplement';t0=time.monotonic();prec=1200;LS=LaurentSeriesRing(GF(2),'q',default_prec=prec);q=LS.gen();s1=LS(0);a6=LS(0)
for n in range(1,prec//83+1):
 z=q**(83*n);s1+=GF(2)(n)*z/(1+z);a6+=GF(2)((5*n**3+7*n**5)//12)*z/(1+z)
def units(aa,x,y):
 a1,a2,a3,a4,a6=aa;aa3=2*y+a3+a1*x;aa2=3*x+a2;aa4=3*x*x+2*x*a2+a4-a1*y;aaa1=2*aa4/aa3+a1;aaa2=aa2-(aa4/aa3)**2-a1*aa4/aa3;b=-aaa2**3/aa3**2;c=-(aaa1*aaa2-aa3)/aa3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return((s-r)/t,t/(s*s-s-r+1))
records=[]
for a in range(1,42):
 uk=q**a;xp=uk/(1+uk)**2;yp=uk**2/(1+uk)**3+s1
 for n in range(1,prec//83+2):
  if 83*n+a<prec:z=q**(83*n+a);xp+=z/(1+z)**2;yp+=z*z/(1+z)**3
  if 83*n-a<prec:w=q**(83*n-a);xp+=w/(1+w)**2;yp+=w/(1+w)**3
 assert(yp*yp+xp*yp+xp**3+s1*xp+a6).valuation()>=prec
 x,y=units([LS(1),LS(0),LS(0),s1,a6],xp,yp);vx,vy=int(x.valuation()),int(y.valuation());assert x.prec()>200 and y.prec()>200;records.append([a,vx,vy])
px=sum(max(0,-x)for a,x,y in records);py=sum(max(0,-y)for a,x,y in records);assert(px,py)==(90,95)
# At the41 width1 geometric cusps, Frobenius acts transitively up to+-.
# ord83(2)=82, so Phi83 is irreducible and zeta is represented exactly.
assert pow(2,41,83)==82 and all(pow(2,k,83)!=1 for k in range(1,82));R=PolynomialRing(GF(2),'z');phi=sum(R.gen()**i for i in range(83));assert phi.is_irreducible();F=GF(2**82,'z',modulus=phi);z=F.gen();assert z**83==1 and z!=1;xp=z/(1+z)**2;yp=z*z/(1+z)**3;xx,yy=units([F(1),F(0),F(0),F(0),F(0)],xp,yp);assert xx and yy and xx**(2**41)==xx and yy**(2**41)==yy
# The rational change between Tate b,c and the stated units is birational.
A=PolynomialRing(QQ,names=('x','y'));K=FractionField(A);x,y=map(K,A.gens());r=(x*x*y-x*y+y-1)/x/(x*y-1);s=(x*y-y+1)/x/y;b=r*s*(r-1);c=s*(r-1);rr=b/c;ss=c*c/(b-c);tt=rr*ss-2*rr+1;assert(ss-rr)/tt==x and tt/(ss*ss-ss-rr+1)==y
out={'status':'PASS','level':83,'Tate_precision':prec,'width83_cusps':41,'width83_orders':records,'width1_cusps':41,'width1_orders':[0,0],'primitive83root_field_degree':82,'width1_residue_degree_mod_sign':41,'pole_degrees':[px,py],'birational_Tate_bc_unit_inverse_verified_over_Q':True,'plane_model_used':False,'canonical_relations_used':False,'elapsed_seconds':time.monotonic()-t0}
p=H/'cusp_orders_verified.json'
if '--readonly'in sys.argv:
 old=json.loads(p.read_text());assert{k:v for k,v in old.items()if k!='elapsed_seconds'}=={k:v for k,v in out.items()if k!='elapsed_seconds'}
else:p.write_text(json.dumps(out,indent=2)+'\n')
print('NONIC83_DIRECT_TATE_CUSP_POLE_BOUND_PASS',json.dumps({k:v for k,v in out.items()if k!='width83_orders'},sort_keys=True))
