from sage.all import *
from pathlib import Path
import re,json,collections,sys
p=Path(__file__).resolve().parent/'data';d=p/'rebound';s=(d/'class_arithmetic33.py').read_text();J=lambda f:json.loads(f.read_text())
def assign(name):
 m=re.search(r'\b'+name+r'\s*=\s*(.*?)(?=\n\w+\s*=|\Z)',s,re.S);assert m,name
 return json.loads(re.sub(r',\s*\]',']',m[1].strip().replace('True','true').replace('False','false')))
rel=assign('relations');assert rel==[[15,0,0],[0,42373650,0]];mods=[15,42373650];bind=J(d/'cusp_binding33.json');qc=assign('rational_cusp_generators');assert qc==bind['global_cusp_classes'];degs=bind['global_cusp_degrees'];assert sorted(degs)==[1]*10+[2]*5+[10]*2 and len({tuple(a)for a in bind['global_cusp_signatures']})==17
A=matrix(ZZ,[v[:2]for v in qc]+[[15,0],[0,42373650]]);assert abs(A.row_module().basis_matrix().det())==3;assert all(v[0]%3==0 for v in qc)
add=lambda a,b:((a[0]+b[0])%15,(a[1]+b[1])%42373650,0)
def effective(maxd,items):
 vv=[set()for _ in range(maxd+1)];ct=[0]*(maxd+1);vv[0]={(0,0,0)};ct[0]=1
 for dd,cc in items:
  for n in range(dd,maxd+1):vv[n]|={add(a,cc)for a in vv[n-dd]};ct[n]+=ct[n-dd]
 return vv,ct
qeff,qcounts=effective(9,list(zip(degs,qc)));assert qcounts[9]==144250 and len(qeff[9])==144250 and qeff[9]=={tuple(a)for a in assign('rational_effective_cusps')}
comp,counts=effective(4,[(dd,cc)for dd,cc in zip(degs,qc)if dd<=4]);assert counts==[1,10,60,270,1005]and[list(map(len,comp))][0]==counts
raw=(d/'noncuspidal_places33.jsonseq').read_text();dec=json.JSONDecoder();places=[]
while raw.strip():raw=raw.lstrip();z,i=dec.raw_decode(raw);places.append(z);raw=raw[i:]
fp=PolynomialRing(ZZ,'Z')(J(p/'frobenius_coefficients33.txt'));assert fp.degree()==42 and fp(1)==635604750
tr=[];closed=[]
for n in range(1,10):
 sn=-n*fp[42-n]-sum(fp[42-j]*tr[n-j-1]for j in range(1,n));tr.append(sn);ct=2**n+1-sn-sum(dd*closed[dd-1]for dd in divisors(n)if dd<n);assert ct%n==0;closed.append(ct//n)
assert closed==[10,5,0,0,4,20,30,25,50]
assert collections.Counter(z['degree']for z in places)=={5:2,6:20,7:30,8:25,9:50};assert len({(z['degree'],tuple(z['field_polynomial']),tuple(map(tuple,z['coordinates'])))for z in places})==127
k=GF(2);R=PolynomialRing(k,'T');P=PolynomialRing(k,2,names=('x','y'));f=P((p.parents[1]/'mdmagma'/'models_X1_n'/'FFFc33.txt').read_text().replace('\n',' ').strip().rstrip(';'));nonsingular=0
for z in places:
 h=R(z['field_polynomial']);assert h.degree()==z['degree']and h.is_irreducible();L=k.extension(h,'a');a=L.gen();x,y=[R(v)(a)for v in z['coordinates']];assert f(x,y)==0
 if any(ff(x,y)!=0 for ff in f.gradient()):nonsingular+=1
 v=z['class'];assert len(v)==3 and v[-1]==0
known=[assign(n)for n in['known_septic_plus_degree2cusp','known_octic_plus_cusp','known_nonic_classes']];assert list(map(len,known))==[600,50,10];ks={tuple(v)for vv in known for v in vv};assert len(ks)==660
alltests=assign('tests');assert len(alltests)==155770;by=collections.defaultdict(collections.Counter);candidates=set()
for degree,v,H,Q in alltests:
 v=tuple(v);assert H==(v[0]%3==0)and Q==(v in qeff[9]);by[degree][v]+=1
 if H and not Q:candidates.add(v)
for dd in range(5,10):
 expected=collections.Counter(add(z['class'],c)for z in places if z['degree']==dd for c in comp[9-dd]);assert by[dd]==expected
assert len(by[0])==146260 and sum(by[0].values())==146260 and {v for v in by[0]if v[0]%3==0}==qeff[9]
assert len(candidates)==670 and ks<=candidates;unknown=candidates-ks;assert len(unknown)==10
short=J(p/'short_cusp_classes33.json');assert short['global_cusp_signatures']==bind['global_cusp_signatures']and {tuple(z['finite_class'])for z in short['representatives']}==unknown
for z in short['representatives']:
 co=z['global_cusp_coefficients'];assert len(co)==17 and sum(a*b for a,b in zip(co,degs))==9
 got=[sum(co[i]*qc[i][j]for i in range(17))%mods[j]for j in range(2)]+[0];assert got==z['finite_class']
 assert sum(max(a,0)*b for a,b in zip(co,degs))==z['positive_degree']and z['positive_degree']-z['negative_degree']==9 and co[:2]==[0,0]
out={'complete':True,'finite_class_tests':155770,'all_cusp_classes':146260,'global_effective_cusp_classes':144250,'known_septic_plus_degree2cusp':600,'known_octic_plus_cusp':50,'known_nonic_classes':10,'unknown_classes':10,'all127_finite_point_coordinates_verified':True,'affine_nonsingular_points':nonsingular,'all_degree9_global_cusp_representatives_verified':True,'no_exclusion_asserted_here':True}
if'--check'in sys.argv:assert J(p/'independent_replay33.json')==out
else:(p/'independent_replay33.json').write_text(json.dumps(out,indent=2)+'\n')
print('NONIC33_INDEPENDENT_CLASS_AND_REPRESENTATIVE_REPLAY_PASS',out,flush=True)
