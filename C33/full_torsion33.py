from sage.all import *
from pathlib import Path
import json,sys
p=Path(__file__).resolve().parent/'data';D=json.loads((p.parents[1]/'common'/'data'/'inputs'/'vanhoeij_nonic_input.json').read_text());z=D['examples'][10];assert z['N']==33 and z['index']==11
R=PolynomialRing(QQ,'T');T=R.gen();poly=lambda a:R([QQ(c)for c in a]);f=poly(z['f']);xp=poly(z['x']);yp=poly(z['y']);assert f.degree()==9 and f.is_irreducible();disc=f.discriminant();assert QQ(disc).valuation(5)==0
k=GF(5);P=PolynomialRing(k,'t');t=P.gen();h=t*t+2*t-1;assert h.is_irreducible()and P(f)%h==0;L=k.extension(h,'a');a=L.gen()
def evalred(v):return sum((L(QQ(c).numerator())/L(QQ(c).denominator()))*a**i for i,c in enumerate(v))
x=evalred(z['x']);y=evalred(z['y']);assert x*y and x*x*y-x
r=(x*x*y-x*y+y-1)/(x*x*y-x);s=(x*y-y+1)/(x*y);b=r*s*(r-1);c=s*(r-1);aa=[1-c,-b,-b,L(0),L(0)];E=EllipticCurve(L,aa);assert E.discriminant()!=0
count=1
for xx in L:
 for yy in L:
  if yy*yy+aa[0]*xx*yy+aa[2]*yy==xx**3+aa[1]*xx*xx:count+=1
assert count==33;Q=E(0,0);assert 33*Q==E(0)and 3*Q!=E(0)and 11*Q!=E(0)
out={'complete':True,'source_example':11,'global_point_degree':9,'global_field_polynomial_irreducible':True,'global_field_discriminant_5adic_valuation':0,'prime':5,'residue_polynomial':[4,2,1],'residue_degree':2,'good_reduction':True,'direct_enumerated_affine_pairs':625,'curve_order':33,'point_order':33,'full_global_torsion_order':33,'reason':'Unramified odd-prime good reduction injects torsion; characteristic-zero point order33 is independently verified in the complete native classification.'}
if'--check'in sys.argv:assert json.loads((p/'full_torsion33.json').read_text())==out
else:(p/'full_torsion33.json').write_text(json.dumps(out,indent=2)+'\n')
print('NONIC33_FULL_TORSION_DIRECT_REDUCTION_PASS order33_good_unramified_F25=true',flush=True)
