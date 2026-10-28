from pathlib import Path
from math import gcd
from collections import Counter
import ast,json,sys
from sage.all import GF,PolynomialRing,EllipticCurve
P=Path(__file__).resolve().parents[1]/'data'/'orbit_section_criterion'

def rd(q):return ast.literal_eval(q.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def key(vs,D):
 out=[]
 for j in range(D):out.append(tuple(int(x) for v in vs for x in v._vector_()));vs=[v**5 for v in vs]
 return min(out)
def norm(E,P,Q):
 x=Q[0]-P[0];y=Q[1];return [x**3/y**2,(E.a2()+3*P[0])*x*x/y**2]
def polfield(pol,d):
 if d==1:return GF(5)
 R=PolynomialRing(GF(5),'z');return GF(5**d,'z',modulus=R(pol))
a=rd(P/'kummer_certified_arithmetic.txt');pol,ky,inert,records,orbits,dims=a;F=polfield(pol,8);z=F.gen();b=sum(F(v)*z**j for j,v in enumerate(ky[:8]));c=sum(F(v)*z**j for j,v in enumerate(ky[8:]));E=EllipticCurve(F,[0,c,0,b*(1-b-c),0]);P2=E(0,0);Q=E(b,b);g=norm(E,P2+18*Q,P2+Q);assert g==[b**(5**6),c**(5**6)]
PP,QQ=P2,Q
for j in range(4):PP,QQ=PP+18*QQ,PP+QQ
assert PP==P2 and QQ==Q
R=PolynomialRing(GF(5),'z');K=GF(25,'z',modulus=R(inert[0]));u=sum(K(v)*K.gen()**j for j,v in enumerate(inert[1]));Y=PolynomialRing(K,'Y').gen();assert (Y**4-u).is_irreducible()
# Independent source-to-quotient invariant and Kummer-unit alignment.
x,y=(4*Q)[0],(4*Q)[1];aa3=2*y;aa2=3*x+c;aa4=3*x*x+2*c*x+b*(1-b-c);m=aa4/aa3;aa1=2*m;aa2-=m*m
B9=-aa2**3/aa3**2;C9=1-aa1*aa2/aa3;scale=aa3/aa2;tv=C9*C9/(B9-C9);wv=4*((9*Q)[0]-(18*Q)[0])/(P2[0]-(P2+18*Q)[0]);rho=((18*Q)[0]-x)/scale**2
r=tv*tv-tv+1;assert C9==tv*(r-1) and B9==r*C9
A2=-B9+(1-C9)**2/4;A4=-(1-C9)*B9/2;beta=A4+2*A2*rho+3*rho*rho;fv=16*beta**3*(1+4/wv)**2
ipol,ivals=rd(P/'inert_basepoint_values.txt');RF=PolynomialRing(F,'v');roots=RF(ipol).roots(multiplicities=False);matched=[]
for root in roots:
 vals=[sum(F(a)*root**j for j,a in enumerate(coeff)) for coeff in ivals]
 if vals[:2]==[tv,wv]:matched.append(vals);assert vals[2:]==[rho,fv]
assert len(matched)==1
assert fv==sum(F(a)*roots[0]**j for j,a in enumerate(inert[1])) or fv==sum(F(a)*roots[1]**j for j,a in enumerate(inert[1]))
source=[];split=[]
for rec in records:
 i,d,v,mod,coeff,dd,hh=rec;assert v%4==0;L=polfield(mod,d);w=L.gen() if d>1 else L(1);u=sum(L(a)*w**j for j,a in enumerate(coeff));RR=PolynomialRing(L,'Y');Y=RR.gen();fs=(Y**4-u).factor();actual=sorted([d*h.degree() for h,m in fs for _ in range(m)]);assert actual==sorted(dd);source+=dd
 if d==1 and actual==[1,1,1,1]:assert u==1 and hh==[1,0,0,0];split.append(i)
assert Counter(source)=={1:24,2:8,6:4} and len(split)==6 and dims==[1]*6 and [r[0] for r in orbits]==split
# Independent full Tate cusp action for g(P,Q)=(P+18Q,P+Q).
n=36;labels={(a,b,c,d) for a in range(n) for b in range(n) if gcd(gcd(a,b),n)==1 for c in range(2) for d in range(2) if (c*b-d*a)%2==1}
def orb(x):
 a,b,c,d=x;return frozenset(((s*(a+k*b))%n,s*b%n,s*(c+k*d)%2,s*d%2) for s in [-1,1] for k in range(n))
cs=sorted(set(map(orb,labels)),key=lambda x:sorted(x));ix={x:i for i,C in enumerate(cs) for x in C};assert len(cs)==64;gp=[];fp=[]
for C in cs:
 gg={ix[((a+18*c)%36,(b+18*d)%36,(c+a)%2,(d+b)%2)] for a,b,c,d in C};ff={ix[(5*a%36,b,c,d)] for a,b,c,d in C};assert len(gg)==len(ff)==1;gp.append(gg.pop());fp.append(ff.pop())
assert all(gp[gp[gp[gp[i]]]]==i for i in range(64))
assert all(gp[i]!=i and gp[gp[i]]!=i and gp[gp[gp[i]]]!=i for i in range(64))
rat={i for i in range(64) if fp[i]==i};assert len(rat)==24;left=set(rat);rr=[]
while left:
 i=min(left);o={i,gp[i],gp[gp[i]],gp[gp[gp[i]]]};assert len(o)==4 and o<=rat;left-=o;rr.append(o)
assert len(rr)==6
out=dict(complete=True,degree=9,excluded_support=[8,1],source_key=ky,source_field_polynomial=pol,source_stabilizer_frobenius_phase=6,cusp_orbits=6,rational_residue_cusps=24,all_cusp_orbit_sizes=4,kummer_inert_degree=4,base_noncuspidal_place_degree=2,source_place_degree=8,independent_source_to_kummer_invariant_alignment=True,source_cusp_signature={'1':24,'2':8,'6':4},rr_character_dimensions=[[1,0,0,0] for _ in range(6)],rr_source_dimensions=[1]*6,necessary_dimension=4)
if '--readonly' in sys.argv:assert json.loads((P/'independent_orbit_certificate.json').read_text())==out
else:(P/'independent_orbit_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C236_INDEPENDENT_G4_TATE_KUMMER_LOCAL_SUPPORT_AND_SECTION_REPLAY_PASS',json.dumps(out),flush=True)
