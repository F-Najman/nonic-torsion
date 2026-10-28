from pathlib import Path
import ast,re,json,time,sys,math
from sage.all import GF,PolynomialRing,PowerSeriesRing,EllipticCurve,Matrix,vector,lcm
P=Path(__file__).resolve().parents[1]/'data'/'CM6';t0=time.monotonic();R0=PolynomialRing(GF(2),'a');fields={};rankcache={};checked={};ranktimes=[];pointcount=0

def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']'))
def rd(f):return lit((P/f).read_text())
def mat(s,name):
 m=re.search(name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,name
 return Matrix(GF(2),int(m[1]),int(m[2]),lit(m[3]))
s=(P.parent/'source_units'/'canonical_unit_data.m').read_text();RX=mat(s,'UnitXRelations');RY=mat(s,'UnitYRelations')

def field(pol):
 ky=tuple(pol)
 if ky not in fields:fields[ky]=GF(2**(len(pol)-1),'a',modulus=R0(pol),impl='pari_ffelt')
 return fields[ky]
def dec(K,c):return sum(K(v)*K.gen()**j for j,v in enumerate(c))
def units(aa,x,y):
 a1,a2,a3,a4,a6=aa;b3=2*y+a3+a1*x;b2=3*x+a2;b4=3*x*x+2*x*a2+a4-a1*y;a1p=2*b4/b3+a1;a2p=b2-(b4/b3)**2-a1*b4/b3;b=-a2p**3/b3**2;c=-(a1p*a2p-b3)/b3;r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)]
def psi(x,aa):
 a1,a2,a3,a4,a6=aa;b2=a1*a1+4*a2;b4=a1*a3+2*a4;b6=a3*a3+4*a6;b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4;s2=a1*x+a3;s3=3*x**4+b2*x**3+3*b4*x*x+3*b6*x+b8;s4=s2*(2*x**6+b2*x**5+5*b4*x**4+10*b6*x**3+10*b8*x*x+(b2*b8-b4*b6)*x+b4*b8-b6*b6);pp=[x.parent()(x.base_ring().zero()),x.parent()(x.base_ring().one()),s2,s3,s4]
 for n in range(5,66):
  m=n//2;v=pp[m+2]*pp[m]**3-pp[m-1]*pp[m+1]**3 if n%2 else pp[m]/s2*(pp[m+2]*pp[m-1]**2-pp[m-2]*pp[m+1]**2);pp.append(v)
 return pp[65]
def pairindex(i,j):
 i,j=sorted([i,j]);return i*85-i*(i-1)//2+j-i

def check_record(rec,label):
 global pointcount
 cachekey=json.dumps(rec,separators=(',',':'))
 if cachekey in checked:return checked[cachekey]
 pol,aa0,p0,z0,xy0,mode,ej,zj=rec;K=field(pol);aa=[dec(K,c)for c in aa0];E=EllipticCurve(K,aa);pp=E(*[dec(K,c)for c in p0]);assert 65*pp==E(0)and 13*pp!=E(0)and 5*pp!=E(0)
 xy=units(aa,*pp.xy());assert xy==[dec(K,v)for v in xy0];degree=int(lcm([x.minpoly().degree()for x in xy]));PS=PowerSeriesRing(K,'t',default_prec=16);t=PS.gen();xx,yy=[PS([dec(K,v)for v in row]).add_bigoh(16)for row in ej];aas=list(map(PS,aa));assert mode==(1 if aa[0]==0 else 6)
 aas[0 if mode==1 else 4]+=t;assert psi(xx,aas).valuation()>=16,(label,psi(xx,aas).valuation(),psi(xx,aas).precision_absolute());assert psi(PS(pp[0])+t,list(map(PS,aa)))[1]!=0,label
 assert (yy*yy+aas[0]*xx*yy+aas[2]*yy-xx**3-aas[1]*xx**2-aas[3]*xx-aas[4]).valuation()>=16 and aa[0]*pp[0]+aa[2]!=0
 U,V=units(aas,xx,yy);zz=[PS([dec(K,v)for v in row]).add_bigoh(16)for row in zj];assert[zz[i][0]for i in range(121)]==[dec(K,v)for v in z0]
 z=zz;assert all(v.valuation()>=0 for v in z);ki=next(i for i,v in enumerate(z)if v[0]);assert (z[ki]-1).valuation()>=16
 # All R1 unit identities plus rank120 identify the genuine canonical arc.
 for M,u in[(RX,U),(RY,V)]:
  assert all(sum(z[j]*(int(v[j])+u*int(v[121+j]))for j in range(121)).valuation()>=16 for v in M.rows())
 M0=Matrix(K,[[K(int(v[j]))+u*int(v[121+j])for j in range(121)]for M,u in[(RX,xy[0]),(RY,xy[1])]for v in M.rows()]);rk=M0.rank();assert rk==120
 assert any(v[1]for v in z)and any(v[1]for v in zz)
 pointcount+=1;print('C65_CM6_GENUINE_R1_JETS_PASS',label,'ambient',K.degree(),'closed',degree,'rank',rk,'seconds',round(time.monotonic()-t0,2),flush=True)
 result=(K,E,pp,degree);checked[cachekey]=result;return result

def el_from(rec):
 K=field(rec[0]);return K,EllipticCurve(K,[dec(K,v)for v in rec[1]]),[dec(K,v)for v in rec[2]]
casecounts=[];orbits=0;geometric=0
for mode in[6]:
 data=rd(f'good{mode}_jets_data.txt');assert data[:2]==[mode,16];_,_,sources,targets,kernels=data
 for ii,(ss,ts,krs)in enumerate(zip(sources,targets,kernels),1):
  K,E,pp,d=check_record(ss[0],f'{mode}/{ii}/P');assert d==mode
  for k,rec in zip([3,9],ss[1:]):
   L,EE,pt,dd=check_record(rec,f'{mode}/{ii}/{k}P');assert L is K and EE==E and pt==k*pp and dd==d
  common=None;norms=[];normdegree=0
  for jj,(trr,kr)in enumerate(zip(ts,krs),1):
   e,td,w,tr,tr5=trr;pol,els,pls,kpc,egs,pgs,ke=kr;assert ke==e and mode*e==td*w;L=field(pol);EL=EllipticCurve(L,[dec(L,v)for v in els]);PL=EL(*[dec(L,v)for v in pls]);EG=EllipticCurve(L,[dec(L,v)for v in egs]);PG=EG(*[dec(L,v)for v in pgs]);RR=PolynomialRing(L,'x');kp=RR([dec(L,v)for v in kpc]);assert kp.degree()==1
   phi=EL.isogeny(kp,degree=3,check=True);assert phi.codomain()==EG and phi(PL)==PG and L.degree()==mode*e
   LL,ET,PT,dd=check_record(tr,f'{mode}/{ii}/T{jj}');LL5,ET5,PT5,dd5=check_record(tr5,f'{mode}/{ii}/5T{jj}');assert LL is LL5 and ET==ET5 and PT5==3*PT and dd==dd5==td
   assert LL is L
   if EG.a1():
    assert EG.a1()==1 and EG.a3()==0;assert ET.a_invariants()==(L(1),EG.a2(),L(0),L(0),EG.a6()+EG.a4()**2)and PT[0]==PG[0]and PT[1]==PG[1]+EG.a4()
   else:assert ET==EG and PT==PG
   if common is None:common=L;EC=EL;PC=PL;mp=lambda v:v
   else:
    roots=PolynomialRing(common,'w')(pol).roots(multiplicities=False);mp=None
    for root in roots:
     fun=lambda v,root=root:sum(common(c)*root**j for j,c in enumerate(v._vector_()))
     if [fun(v)for v in EL.a_invariants()]==list(EC.a_invariants())and[fun(v)for v in PL.xy()]==list(PC.xy()):mp=fun;break
    assert mp is not None
   RC=PolynomialRing(common,'x');gg=RC([mp(v)for v in kp]);assert next(a for a in range(1,5)if all(c**(2**(mode*a))==c for c in gg))==e
   ng=RC(1)
   for a in range(e):ng*=RC([c**(2**(mode*a))for c in gg])
   norms.append(ng);orbits+=1;geometric+=e;normdegree+=td*w
  assert sum(tr[0]for tr in ts)==4 and normdegree==4*mode
  assert all(norms[i].gcd(norms[j]).degree()==0 for i in range(len(norms))for j in range(i))
  prod=PolynomialRing(common,'x')(1)
  for ng in norms:prod*=ng
  expected=PolynomialRing(common,'x')(EC.division_polynomial(3));expected/=expected.leading_coefficient();assert prod==expected
  casecounts.append(dict(source_degree=mode,source_index=ii,Hecke_closed_degrees=[r[1]for r in ts],Hecke_norm_weights=[r[2]for r in ts]));print('C65_CM6_ALL_FOUR_KERNELS_AND_LEVEL_IMAGES_PASS',mode,ii,flush=True)
out=dict(complete=True,level=65,prime=2,genus=121,canonical_jet_precision=16,independent_genuine_R1_arcs=pointcount,source_cases=casecounts,independent_kernel_orbits=orbits,geometric_cyclic3_kernels=geometric,all_division_polynomial_norms_exact=True,all_actual_P_3P_9P_images_verified=True,supersingular_a1_deformation=True,all_unit_recovery_ranks120=True,all_canonical_tangents_nonzero=True)
if '--readonly'in sys.argv:assert json.loads((P/'cm6_sage_certificate.json').read_text())==out
else:(P/'cm6_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C65_CM6_GOOD_HECKE_AND_GENUINE_JETS_INDEPENDENT_SAGE_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
