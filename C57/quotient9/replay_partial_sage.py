from pathlib import Path
import ast,re,json,sys,time
from collections import Counter,deque
from sage.all import GF,PolynomialRing,EllipticCurve_from_j,Matrix,divisors
ROOT=Path(__file__).resolve().parents[2];P=ROOT/'C57'/'data'/'quotient9';t0=time.monotonic()
def literal(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
def rd(f):return literal((P/f).read_text())
def norm(z):return [v/next(w for w in z if w)for v in z]
def fkey(z,d):return min(tuple(int(c)for v in norm([v**(2**a)for v in z])for c in v._vector_())for a in range(d))
R=PolynomialRing(GF(2),'x')
canonical=(P/'canonical_data.m').read_text();d5=(P/'diamond5_data.m').read_text()
def mat(s,name):
 m=re.search(name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,name
 return Matrix(GF(2),int(m[1]),int(m[2]),literal(m[3]))
Q=mat(canonical,'qcoeffs');rel=mat(canonical,'QuadRelations');D=mat(d5,'Diamond5Matrix');assert D*D==Matrix.identity(GF(2),9)and(D-Matrix.identity(GF(2),9)).rank()==4
pairs=[(i,j)for i in range(9)for j in range(i,9)]
PS=__import__('sage.all',fromlist=['PowerSeriesRing']).PowerSeriesRing(GF(2),'q',default_prec=120);q=PS.gen();forms=[PS(list(v)[:120]).add_bigoh(120)for v in Q.rows()];products=[forms[i]*forms[j]for i,j in pairs]
assert all(sum(c*v for c,v in zip(row,products))==0 for row in rel.rows())
full,kept,allplaces=rd('noncuspidal_filter_v2_data.txt');bydegree=Counter();points={};seen={}
for e,i,pol,co,cusp in allplaces:
 K=GF(2)if e==1 else GF(2**e,'a',modulus=R(pol));a=K.gen();z=[sum(K(c)*a**k for k,c in enumerate(v))for v in co];assert z==norm(z)
 assert min(d for d in divisors(e)if norm([v**(2**int(d))for v in z])==z)==e
 prod=[z[i]*z[j]for i,j in pairs];assert all(sum(K(c)*v for c,v in zip(row,prod))==0 for row in rel.rows())
 ky=fkey(z,e);assert(e,ky)not in seen;seen[e,ky]=i;points[e,i]=z;bydegree[e]+=1
fp=rd('quotient_frobenius2.txt');assert len(fp)==19 and fp[-1]==1 and fp[0]==2**9;tr=[];nc=[]
for n in range(1,10):
 s=-n*fp[18-n]-sum(fp[18-j]*tr[n-j-1]for j in range(1,n));tr.append(s);v=2**n+1-s-sum(int(d)*nc[int(d)-1]for d in divisors(n)if d<n);assert v%n==0;nc.append(v//n)
assert [bydegree[e]for e in range(1,10)]==nc==[6,9,0,2,0,0,24,24,64]
# Independent complete enumeration of j and twist orders on the source fields.
jsets={};curvecounts={}
for e in range(6,10):
 txt=(ROOT/f'C57/data/census/N57_f{e}_records.m').read_text();pol=literal(re.search(r'FieldPolynomial:=(\[.*?\]);',txt,re.S)[1]);expected={tuple(v)for v in literal(re.search(r'QualifiedJPolynomials:=(\[.*?\]);',txt,re.S)[1])};K=GF(2**e,'a',modulus=R(pol));js=set();qualified=0;tested=0
 for j in K:
  for E in EllipticCurve_from_j(j).twists():
   tested+=1
   if E.cardinality()%57==0:qualified+=1;js.add(tuple(int(c)for c in j.minpoly()))
 assert js==expected and qualified==int(re.search(r'QualifiedCurves:=(\d+);',txt)[1]);jsets[e]=js;curvecounts[e]=dict(twists=tested,qualified=qualified,j_polynomials=sorted(map(list,js)))
 print('C57_COMPLETE_SOURCE_TWIST_CENSUS_PASS',e,tested,qualified,flush=True)
# Reconstruct the finite group, alpha image and exact dual-A5 kernel by BFS.
g,rels,cd,widths,sig,sf,hc,ac,jci,ai,kg,a5=rd('finite_group_A5_dual_verified_data.txt');assert g==[105,210,0]and ai==[7,7]
def key(v):assert len(v)==3 and v[-1]==0;return(v[0]%105,v[1]%210)
def add(a,b):return((a[0]+b[0])%105,(a[1]+b[1])%210)
def mul(n,a):return(n*a[0]%105,n*a[1]%210)
ac=list(map(key,ac));a5=list(map(key,a5));hc=list(map(key,hc));operator={(0,0):(0,0)};todo=deque([(0,0)])
while todo:
 v=todo.popleft()
 for aa,bb in zip(ac,a5):
  vv=add(v,aa);ww=add(operator[v],bb)
  if vv in operator:assert operator[vv]==ww
  else:operator[vv]=ww;todo.append(vv)
assert len(operator)==49;kernel={v for v,w in operator.items()if w==(0,0)};assert len(kernel)==7
ks={(0,0)}
for gen in map(key,kg):ks={add(v,mul(n,gen))for v in ks for n in range(7)}
assert ks==kernel
assert all(mul(9,ac[i])not in kernel for i in [2,4])
cusp=json.loads((P/'cusp57.json').read_text());assert cusp['all_source_cusps_degree_le3_globally_intact'] and cusp['only_non_global_source_cusps_degree_le9_two_single_nonics']
assert widths==[1,19,57,3,57,1]and sig==[6,2,5,4,3,1]and sf==[1,2,5,4,3,6]
actual=[]
for row in full:
 d,e,w,jp,i,si,cl,yes,pol,co=row;assert d%e==0 and w==d//e and tuple(jp)in jsets[d]
 z=points[e,i];sz=norm([sum(int(D[a,b])*z[b]for b in range(9))for a in range(9)]);assert fkey(sz,e)==fkey(points[e,si],e)
 assert bool(yes)==(mul(w,key(cl))in kernel)
 if yes:actual.append(row)
assert actual==kept and len(full)==23 and[(r[0],r[1],r[4],r[5])for r in kept]==[(6,2,5,5),(8,8,8,19),(8,8,19,8)]
out=dict(complete=True,partial_only=True,degree=9,level=57,quotient_genus=9,prime=2,independent_canonical_points=len(allplaces),closed_place_counts=nc,source_census=curvecounts,finite_alpha_order=49,dual_A5_kernel_order=7,tested_conservative_non_cusp_cases=23,surviving_images=[dict(source_degree=r[0],target_degree=r[1],norm=r[2],target_place=r[4],diamond5_place=r[5])for r in kept],all_source_degree7_and_degree9_noncuspidal_reductions_excluded=True,non_global_entirely_cuspidal_reductions_excluded=True,entirely_cuspidal_reductions_open=True,full_C57_exclusion=False)
if '--readonly'in sys.argv:assert json.loads((P/'partial_sage_certificate.json').read_text())==json.loads(json.dumps(out))
else:(P/'partial_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C57_QUOTIENT9_PARTIAL_INDEPENDENT_SAGE_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
