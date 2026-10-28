from pathlib import Path
import ast,json,time,hashlib,sys
from sage.all import GF,PolynomialRing,Matrix,vector,lcm
P=Path(__file__).resolve().parents[1]/'data'/'CM6';t0=time.monotonic();R=PolynomialRing(GF(2),'a');keylist=[];places=[];ids={};cases=[];cache={}
def rd(name):return ast.literal_eval((P/name).read_text().replace('[*','[').replace('*]',']'))
def point(rec,ref):
 raw=json.dumps(rec[0:5],separators=(',',':'))
 if raw in cache:return cache[raw]
 pol,aa,pp,co,xy=rec[:5];m=len(pol)-1;K=GF(2**m,'a',modulus=R(pol));a=K.gen();z=[sum(K(c)*a**j for j,c in enumerate(v))for v in co];pivot=next(v for v in z if v);z=[v/pivot for v in z];d=int(lcm([v.minpoly().degree()for v in z]));kk=next(i for i,v in enumerate(z)if v.minpoly().degree()==d);u=z[kk];pm=list(map(int,u.minpoly()));B=Matrix(GF(2),[list((u**j)._vector_())for j in range(d)]);C=Matrix(GF(2),[list(v._vector_())for v in z]);cp=B.solve_left(C);assert cp*B==C;key=[pm,[list(map(int,v))for v in cp.rows()]];s=json.dumps(key,separators=(',',':'));h=hashlib.sha256(s.encode()).hexdigest()
 if s not in ids:
  ids[s]=len(places)+1;keylist.append(dict(id=ids[s],sha256=h,closed_degree=d,primitive_coordinate=kk+1,key=key));places.append([d,m,rec[7]])
 ix=ids[s];cache[raw]=(ix,d);return ix,d
for mode in[6]:
 data=rd(f'good{mode}_jets_data.txt');assert data[:2]==[mode,16];_,_,sources,targets,kernels=data
 for variant,weights in enumerate([[1,2,-3],[3,-2,-1]],1):
  ii=variant;ss=sources[0];ts=targets[0]
  terms={}
  def add(rec,n,ref):
   i,d=point(rec,ref);terms[i]=terms.get(i,0)+n;return d
  for j,n in enumerate(weights):assert add(ss[j],n,['S',mode,ii,j+1])==mode
  for j,(e,d,w,tr,tr5)in enumerate(ts,1):
   assert mode*e==d*w and add(tr,-w,['T',mode,ii,j,1])==add(tr5,w,['T',mode,ii,j,5])==d
  terms={i:n for i,n in terms.items()if n};pos=sum(places[i-1][0]*n for i,n in terms.items()if n>0);neg=sum(places[i-1][0]*(-n)for i,n in terms.items()if n<0);assert pos==neg
  assert all(abs(n)<=16 for n in terms.values());cases.append([mode,ii,[[i,n]for i,n in sorted(terms.items())]])
  print('C65_NECESSARY_RR_DIVISOR',mode,ii,'positive_degree',pos,'signed_multiplicities',list(terms.values()),flush=True)
out=dict(complete=True,prime=2,genus=121,line_bundle='K',definition='point A3=T3-1-3h; variant1 necessary, variant2 extra',closed_place_keys=keylist,case_count=len(cases),unique_places=len(places),cases=cases)
if '--readonly' in sys.argv:assert rd('rr_places_data.txt')==places and rd('rr_divisors_data.txt')==cases and json.loads((P/'rr_divisors_sage_certificate.json').read_text())==out
else:(P/'rr_places_data.txt').write_text(str(places).replace('[','[*').replace(']','*]')+'\n');(P/'rr_divisors_data.txt').write_text(str(cases).replace('[','[*').replace(']','*]')+'\n');(P/'rr_divisors_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C65_NECESSARY_RR_DIVISORS_ASSEMBLED',len(cases),len(places),'seconds',time.monotonic()-t0,flush=True)
