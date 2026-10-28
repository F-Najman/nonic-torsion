from pathlib import Path
import ast,json,re,sys,time
from collections import Counter
from sage.all import GF,QQ,ZZ,PolynomialRing,EllipticCurve,divisors,matrix,vector
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
S=Path(__file__).resolve().parents[1]/'data';P=S/'genus5_quotient';t0=time.monotonic()
def rd(p):
 s=p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False')
 s=re.sub(r'\([^()\n]*F\.1[^()\n]*\)',lambda m:json.dumps(m.group()),s)
 return ast.literal_eval(s)
class Quotient:
 def __init__(self,R,H):
  self.n=len(R[0])-1;self.M=hermite_normal_form(Matrix([r[:-1] for r in R+H]).T);self.index=abs(int(self.M.det()))
 def key(self,v):
  assert len(v)==self.n+1 and v[-1]==0;w=list(v[:-1])
  for j in range(self.n-1,-1,-1):
   q=w[j]//int(self.M[j,j])
   for i in range(j+1):w[i]-=q*int(self.M[i,j])
  return tuple(w)
def add(*vs):return [sum(x)for x in zip(*vs)]
def mul(n,v):return [n*x for x in v]
def projection(E,P2,Q):
 S7=6*Q;T3=14*Q;R2=21*Q;Z2=P2+R2;x,y=S7[0],S7[1];m=(3*x*x+2*E.a2()*x+E.a4())/(2*y);a2p=3*x+E.a2()-m*m;a3p=2*y
 k=a3p/a2p;b=-a2p**3/a3p**2;c=1-2*m*a2p/a3p;v=b/c;assert b==v*v*(v-1)and c==v*(v-1);H=v*v-v+1;assert H
 u=4*(P2[0]-R2[0])*(P2[0]-Z2[0])*(R2[0]-Z2[0])/(k**6*v**4*(v-1)**4);x3=(T3[0]+E.a2()/3)/(k*k*H*H);y3=T3[1]/(k**3*H**3)
 t=u*u;D=t*t+13*t+49;N=t*t+5*t+1;M=t**4+14*t**3+63*t*t+70*t-7
 assert y3*y3==x3**3-N/(48*D)*x3+M/(864*D*D)
 return [u,x3,y3]
def degree(vs,d):return min(int(e)for e in divisors(d)if all(v**(5**int(e))==v for v in vs))
R5=PolynomialRing(GF(5),'z');fcache={}
def fld(pol):
 k=tuple(pol)
 if k not in fcache:fcache[k]=GF(5**(len(pol)-1),'z',modulus=R5(pol))
 return fcache[k]
def ev(K,cs):return sum(K(c)*K.gen()**i for i,c in enumerate(cs))
def curve(d):
 pol,ai,pv,qv=d;K=fld(pol);E=EllipticCurve(K,[ev(K,a)for a in ai]);pp=[ev(K,a)for a in pv];qq=[ev(K,a)for a in qv];assert pp[2]==qq[2]==1;return K,E,E(pp),E(qq)
# Rational descent identities, checked over Q(v), independent of both native finite probes.
RV=PolynomialRing(QQ,'v');V=RV.fraction_field();v=V.gen();b=v*v*(v-1);c=v*(v-1);H=v*v-v+1;s=(v**3-3*v+1)/(v*(v-1));t=s-8;D=t*t+13*t+49;N=t*t+5*t+1;M=t**4+14*t**3+63*t*t+70*t-7
ET=EllipticCurve(V,[1-c,-b,-b,0,0]);assert D==H**3/(v*v*(v-1)**2)and ET.c4()/H**4==N/D and ET.c6()/H**6==-M/(D*D)and ET.discriminant()/H**12==t/D**4
assert ET.discriminant()==t*v**8*(v-1)**8
assert ((v-1)/v)**2-((v-1)/v)+1==H/v**2
print('C242_G5_Q_DESCENT_IDENTITIES_SAGE_PASS',flush=True)
# Link the complete Hecke export to exactly the same native group presentation as the filter.
ha=rd(P/'audit_hecke/arithmetic.txt');hf=rd(P/'filter_919/arithmetic.txt');assert ha==hf
iv,rels,hc,acs,records,kgens,kinv,aim=ha;assert iv==[2,6,6,48,0]and kinv==[2,2,2,24]and aim==[3,6]
gp=[r[1]for r in records];ap=[r[3]for r in records];G=Quotient(rels,[]);assert G.index==3456
Smat=matrix(ZZ,[v[:-1]for v in hc+gp]).transpose();Tmat=matrix(ZZ,[v[:-1]for v in acs+ap]).transpose();Rmat=matrix(ZZ,[v[:-1]for v in rels]).transpose();assert Quotient(rels,hc+gp).index==1
for v in Smat.augment(-Rmat).right_kernel().basis():assert G.key(list(Tmat*vector(ZZ,v[:len(hc)+len(gp)]))+[0])==(0,0,0,0)
kcalc=[list(Smat*vector(ZZ,v[:len(hc)+len(gp)]))+[0]for v in Tmat.augment(-Rmat).right_kernel().basis()]
# Reconstruct the exact ambient kernel from the full relation-compatible operator.
# The historical native subgroup-generator dump used relative subgroup coordinates.
KQ=Quotient(rels,kcalc);assert KQ.index==18
KL=matrix(ZZ,[[int(a)for a in r]for r in KQ.M.tolist()]);KI=matrix(ZZ,KL.inverse()*Rmat);assert [abs(int(a))for a in KI.smith_form()[0].diagonal()if abs(a)>1]==kinv
# Native cusp ordering p5: torus at0,torus at0,node0,nodeinfty-torus,nodeinfty,nodeinfty.
p5=rd(P/'p5_arithmetic.txt');branches={i:[j for j,r in enumerate(p5[4])if r[4]==i]for i in range(4)};assert p5[7]==hc
sigma=list(range(len(hc)));delta=list(range(len(hc)));a,b=branches[0];sigma[a],sigma[b]=b,a
for bb in [0,2]:
 a,b=branches[bb];delta[a],delta[b]=b,a
for i in range(len(hc)):
 expected=add(hc[sigma[i]],mul(11,hc[delta[sigma[i]]]),mul(-1,hc[i]),mul(-11,hc[delta[i]]));assert G.key(acs[i])==G.key(expected)
assert all(KQ.key(v)==(0,0,0,0)for v in p5[9])
# Independently reconstruct every cyclic11 kernel isogeny with Sage's checked polynomial isogeny implementation.
proof=rd(P/'audit_hecke/explicit_hecke_data.txt');assert len(proof)==3;isogcount=0
for i,(src,vals,imgs)in enumerate(proof):
 K,E,P2,Q=curve(src);d=K.degree();assert d==9 and Q.order()==42 and P2.order()==2 and P2!=21*Q
 pv=projection(E,P2,Q);assert degree(pv,d)==3
 assert projection(E,P2,11*Q)==[pv[0],pv[1],-pv[2]]
 RK=PolynomialRing(K,'w');normprod=RK.one();multiplicities=Counter()
 for e,srcL,gdat,outL,embed,expectedvals in imgs:
  L,EL,PL,QL=curve(srcL);LL,EG,PG,QG=curve(outL);assert L is LL;z=ev(L,embed);assert R5(src[0])(z)==0;em=K.hom([z],L)
  assert [em(a)for a in E.a_invariants()]==list(EL.a_invariants())and EL([em(a)for a in P2])==PL and EL([em(a)for a in Q])==QL
  RL=PolynomialRing(L,'w');g=RL([ev(L,a)for a in gdat]);assert g.degree()==5 and g.is_monic()
  iso=EL.isogeny(g,check=True);assert iso.degree()==11 and iso.codomain()==EG and iso(PL)==PG and iso(QL)==QG
  vv=projection(EG,PG,QG);assert vv==[ev(L,a)for a in expectedvals]
  ng=RL.one();gg=g
  for j in range(e):ng*=gg;gg=RL([a**(5**d)for a in gg])
  assert gg==g;B=matrix(GF(5),[(z**j)._vector_()for j in range(d)]);cs=[]
  for a in ng:
   co=B.solve_left(vector(GF(5),a._vector_()));b=ev(K,list(co));assert em(b)==a;cs.append(b)
  normprod*=RK(cs);isogcount+=1
  md=degree(vv,L.degree());assert d*e%md==0
  keys=[];vv0=vv
  for j in range(L.degree()):keys.append(tuple(int(a)for v in vv0 for a in v._vector_()));vv0=[v**5 for v in vv0]
  multiplicities[(md,min(keys))]+=d*e//md
 psi=RK(E.division_polynomial(11)).monic();assert normprod==psi and sum(e for e,*_ in imgs)==12
 assert all(m%3==0 for m in multiplicities.values())and sum(dg*m//3 for(dg,_),m in multiplicities.items())==36
 print('C242_G5_HECKE11_SAGE_GENERATOR_PASS',i+1,'isogenies',len(imgs),flush=True)
# Every source projection and finite-field degree in the 919-support filter.
f=rd(P/'filter_919/filter_arithmetic.txt');fiv,frel,kgen,ki,cd,cc,comps,images,tests,left=f;assert fiv==iv and frel==rels and ki==kinv and cc==hc  # old filter stored relative subgroup coordinates; its exact ambient embedding is certified above.
raw=rd(S/'filter_quotient14_linear_arithmetic.txt');own={(r[0],r[1]):r for r in raw[2]};small={(r[0],r[1]):r for r in raw[3]};fields={d:fld(pol)for d,pol in raw[1]};inputleft=rd(S/'filter_quotient14_isogeny_arithmetic.txt')[3];assert len(inputleft)==919
assert len(comps)==6;Cp=[{KQ.key([0]*5)}]+[set()for i in range(5)]
for j,w in enumerate(cd):
 for d in range(w,6):
  for vv in list(Cp[d-w]):Cp[d].add(KQ.key(add(list(vv)+[0],cc[j])))
assert all(len(Cp[i])==len(comps[i])for i in range(6))
cache={};unsupported=0
for ky,supported,e,w,cl in images:
 d=ky[0];K=fields[d];b=ev(K,ky[1:1+d]);c=ev(K,ky[1+d:]);E=EllipticCurve(K,[0,c,0,b*(1-b-c),0]);P2=E(0,0);Q=E(b,b);assert Q.order()==42 and P2!=21*Q
 if not supported:assert E.j_invariant()==0 and e==w==0;unsupported+=1
 else:vv=projection(E,P2,Q);assert degree(vv,d)==e and w==d//e
 cache[tuple(ky)]=(supported,cl)
for j,desc,cdeg,supported,cl,yes in tests:
 assert desc==inputleft[j-1];rows=[own[desc[1],desc[2]]]
 if desc[0]==2:rows.append(small[desc[3],desc[4]])
 assert cdeg==9-sum(r[0]for r in rows)
 bits=[cache[tuple([r[0]]+r[2])]for r in rows];actual=add(*(r[1]for r in bits));assert supported==all(r[0]for r in bits)and G.key(cl)==G.key(actual)
 assert yes==((not supported)or KQ.key(mul(-1,cl))in Cp[cdeg])
assert left==[r[1]for r in tests if r[-1]]and len(left)==119
allowed={tuple(d)for d in json.loads((S/'stabilizer_filter_sage.json').read_text())['remaining']};tail=[d for d in left if tuple(d)in allowed];assert len(tail)==63
out=dict(complete=True,rational_descent_over_Q=True,prime=5,genus=5,finite_jacobian_order=3456,hecke11_kernel_order=192,hecke11_kernel_invariants=kinv,independent_checked_11_isogenies=isogcount,input_supports=919,output_supports=119,after_optional_even_stabilizer=63,source_projections=len(images),unsupported_CM_images_retained=unsupported)
if '--readonly'in sys.argv:assert json.loads((P/'genus5_sage_certificate.json').read_text())==out
else:(P/'genus5_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n');(P/'after_genus5_and_stabilizer.txt').write_text(str(tail)+'\n')
print('C242_GENUS5_INDEPENDENT_SAGE_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
