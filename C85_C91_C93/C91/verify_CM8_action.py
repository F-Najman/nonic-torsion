#!/usr/bin/env python3
"""Exact C91 supersingular degree-eight cusp-complement obstruction."""
from pathlib import Path
import sys,json,time,math,collections
from sage.all import GF,PolynomialRing,EllipticCurve
from unit_arithmetic import *
H=Path(__file__).resolve().parents[1]/'data/C91';R=H.parents[2];s=(R/'C85_C91_C93/data/N91_f8/records.m').read_text();records=lit(s,'Records');assert len(records)==2;pf=lit(s,'FieldPolynomial');PR=PolynomialRing(GF(2),'a');F=GF(256,'a',modulus=PR(pf));FB=Field(bits(pf));enc=lambda a:bits([int(c)for c in a.polynomial()]);dec=lambda a:F(PR(a));reps=[u for u in range(1,46)if math.gcd(u,91)==1]
def units(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P[0],P[1];aa3=a3+a1*x;aa2=x+a2;aa4=x*x+a4+a1*y;w=aa4/aa3;aaa2=aa2+w*w+a1*w;b=aaa2**3/aa3**2;c=(a1*aaa2+aa3)/aa3;r=b/c;ss=c*c/(b+c);v=r*ss+1;xx,yy=(ss+r)/v,v/(ss*ss+ss+r+1);assert xx and yy and xx*yy!=1;ri=(xx*xx*yy+xx*yy+yy+1)/(xx*(xx*yy+1));si=(xx*yy+yy+1)/(xx*yy);assert ri*si*(ri+1)==b and si*(ri+1)==c;return[xx,yy]
allkeys=set();tables=[]
for rec in records:
 E=EllipticCurve(F,list(map(dec,rec[1])));P=E(list(map(dec,rec[2])));assert E.j_invariant()==0 and E.cardinality()==273 and P.order()==91
 def key(P):
  xy=units(E,P);degrees=[d for d in[1,2,4,8]if all(FB.frob(enc(x),d)==enc(x)for x in xy)];assert min(degrees)==8;return FB.key([1]+list(map(enc,xy)),8)
 kk=[key(u*P)for u in reps];assert len(set(kk))==3;assert all(key((27*u)*P)==k for u,k in zip(reps,kk));allkeys.update(kk);tables.append({'curve_coefficients':rec[1],'diamond_representatives':reps,'closed_keys':[[k[0],list(k[1])]for k in kk],'h27_fixes_every_closed_key':True})
assert len(allkeys)==3
out={'status':'PASS','N':91,'CM_j':0,'closed_degree':8,'closed_points':3,'supersingular_census_twists_with91':2,'full_diamond_projective_order':36,'h_scalar':27,'h_order':2,'h_fixes_all_CM8_closed_points':True,'Tate_normal_form_inverse_regular_and_exact_on_all_retained_points':True,'tables':tables};p=H/'CM8_action_verified.json'
if '--readonly'in sys.argv:assert json.loads(p.read_text())==out
else:p.write_text(json.dumps(out,indent=2)+'\n')
print('NONIC91_CM8_ACTION_PASS closed_points3 h27_fixes_all=true')
