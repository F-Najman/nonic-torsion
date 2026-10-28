from pathlib import Path
import time,json,sys
from math import gcd
from sage.all import GF,LaurentSeriesRing,PolynomialRing,QQ,FractionField
P=Path(__file__).resolve().parents[1]/'data'/'source85';t0=time.monotonic();prec=700
def units(aa,x,y):
 a1,a2,a3,a4,a6=aa;aa3=2*y+a3+a1*x;aa2=3*x+a2;aa4=3*x*x+2*x*a2+a4-a1*y;aaa1=2*aa4/aa3+a1;aaa2=aa2-(aa4/aa3)**2-a1*aa4/aa3;b=-aaa2**3/aa3**2;c=-(aaa1*aaa2-aa3)/aa3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return((s-r)/t,t/(s*s-s-r+1))
records=[]
for w in [57,19]:
 F=GF(2)if w==57 else GF(4,'a');LS=LaurentSeriesRing(F,'q',default_prec=prec);q=LS.gen();zeta=F(1)if w==57 else F.gen();assert zeta**(57//w)==1;s1=LS(0);a6=LS(0)
 for n in range(1,prec//w+1):
  z=q**(w*n);s1+=F(n)*z/(1+z);a6+=F((5*n**3+7*n**5)//12)*z/(1+z)
 exponents=[a for a in range(1,29)if gcd(a,57)==1]if w==57 else list(range(1,19))
 for a in exponents:
  u=zeta*q**a;xp=u/(1+u)**2;yp=u*u/(1+u)**3+s1
  for n in range(1,prec//w+2):
   if w*n+a<prec:z=zeta*q**(w*n+a);xp+=z/(1+z)**2;yp+=z*z/(1+z)**3
   if w*n-a<prec:z=zeta**-1*q**(w*n-a);xp+=z/(1+z)**2;yp+=z/(1+z)**3
  assert(yp*yp+xp*yp+xp**3+s1*xp+a6).valuation()>=prec
  x,y=units([LS(1),LS(0),LS(0),s1,a6],xp,yp);assert min(x.prec(),y.prec())>200;records.append(dict(width=w,exponent=a,vx=int(x.valuation()),vy=int(y.valuation())))
 print('C57_SOURCE85_POLE_CUSP_FAMILY_PASS',w,len(exponents),flush=True)
# Width3 and1 cusps have a primitive19- or57-root on the nodal torus.
# Their constant unit values suffice, since both are nonzero.
F=GF(2**18,'z');z=F.multiplicative_generator();constants=[]
for d in [19,57]:
 zz=z**((2**18-1)//d);assert zz.multiplicative_order()==d;aa=[a for a in range(1,(d-1)//2+1)if gcd(a,d)==1]
 for a in aa:
  u=zz**a;xp=u/(1+u)**2;yp=u*u/(1+u)**3;x,y=units([F(1),F(0),F(0),F(0),F(0)],xp,yp);assert x and y
 constants.append(dict(toric_order=d,representatives=len(aa),geometric_cusps=18,unit_valuations=[0,0]))
px=sum(max(0,-r['vx'])for r in records);py=sum(max(0,-r['vy'])for r in records)
assert len(records)==36 and(px,py)==(37,39)
A=PolynomialRing(QQ,names=('x','y'));K=FractionField(A);x,y=map(K,A.gens());r=(x*x*y-x*y+y-1)/x/(x*y-1);s=(x*y-y+1)/x/y;b=r*s*(r-1);c=s*(r-1);rr=b/c;ss=c*c/(b-c);tt=rr*ss-2*rr+1;assert(ss-rr)/tt==x and tt/(ss*ss-ss-rr+1)==y
out=dict(complete=True,level=57,prime=2,Tate_precision=prec,geometric_cusps=72,Tate_nonconstant_cusps=records,toric_constant_cusps=constants,pole_degrees=[px,py],canonical_section_degree_bounds=[168+px,168+py],plane_model_used=False,canonical_relations_used=False,birational_Tate_inverse_over_Q=True)
if '--readonly'in sys.argv:assert json.loads((P/'pole_degrees_sage_certificate.json').read_text())==out
else:(P/'pole_degrees_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C57_SOURCE85_DIRECT_POLE_DEGREES_SAGE_PASS',px,py,'seconds',time.monotonic()-t0,flush=True)
