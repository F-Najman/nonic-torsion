from pathlib import Path
import ast,json,sys
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
P=Path(__file__).resolve().parents[1]/'data'/'diagonal_quotient'

def rd(n):return ast.literal_eval((P/n).read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def fk(vs,D):
 out=[]
 for i in range(D):out.append(tuple(int(x) for v in vs for x in v._vector_()));vs=[v**5 for v in vs]
 return min(out)
pol,ai,pv,qv,key=rd('proof_generator_point_data.txt');D=len(pol)-1;assert D==9;R=PolynomialRing(GF(5),'z');K=GF(5**D,'z',modulus=R(pol));z=K.gen()
def ev(v):return sum(K(a)*z**j for j,a in enumerate(v))
E=EllipticCurve(K,[ev(v) for v in ai]);pp=[ev(v) for v in pv];qq=[ev(v) for v in qv];assert len(pp)==len(qq)==3 and pp[2]==qq[2]==1;P2=E(pp);Q=E(qq)
assert P2==E(0,0) and 2*P2==E(0) and 36*Q==E(0) and 18*Q!=E(0) and 12*Q!=E(0) and P2!=18*Q
Q12=3*Q;T3=12*Q;R9=4*Q;lam=Q12[1]/Q12[0];xp=R9[0]+(R9+T3)[0]+(R9-T3)[0]-2*T3[0];vals=[Q12[0]/lam**2,E.a2()/lam**2,xp/lam**2];assert fk(vals,D)==tuple(key)
def mindeg(vs):return min(int(d) for d in divisors(D) if all(v**(5**int(d))==v for v in vs))
assert mindeg(vals)==3
l36=Q[1]/Q[0];source=[Q[0]/l36**2,E.a2()/l36**2];assert mindeg(source)==9
out=dict(complete=True,prime=5,source_closed_place_degree=9,quotient_closed_place_degree=3,primitive_level=36,full_two_torsion=True,independent_velu3_projection=True,source_field_polynomial=pol)
if '--readonly' in sys.argv:assert json.loads((P/'proof_generator_sage_certificate.json').read_text())==out
else:(P/'proof_generator_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C236_SAGE_EXPLICIT_A7_GENERATOR_SOURCE_AND_PROJECTION_PASS',json.dumps(out),flush=True)
