from pathlib import Path
import ast,re,json,time,sys
from sage.all import GF,Matrix,vector,PolynomialRing,PowerSeriesRing
P=Path(__file__).resolve().parents[1]/'data'/'source85';t0=time.monotonic()
def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']'))
def mat(s,name):
 m=re.search(name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,name
 return Matrix(GF(2),int(m[1]),int(m[2]),lit(m[3]))
s=(P/'qbasis_data.m').read_text();Q=mat(s,'qcoeffs');D=mat(s,'Diamond5Matrix');s=(P/'annihilator_data.m').read_text();W=mat(s,'AnnihilatingDifferentials');assert D**18==Matrix.identity(GF(2),85) and (W*D).row_space()==W.row_space()
ns=lit(re.search(r'HeckeIndices:=(\[.*?\]);',s,re.S)[1]);aq=lit(re.search(r'PositiveFactorEigenvalues:=(\[.*?\]);',s,re.S)[1]);raw=lit(re.search(r'HeckeMatrices:=(\[.*?\]);AnnihilatorMetadata',s,re.S)[1]);T={n:Matrix(GF(2),85,85,r)for n,r in zip(ns,raw)};I=Matrix.identity(GF(2),85);ops=[]
for n,a in zip(ns,aq):
 for ell in[5,7,11,13]:
  kk=[k for k in range(18)if pow(5,k,57)in[ell%57,-ell%57]];assert len(kk)==1;ops+=list(((T[n]-GF(2)(a)*I)*(T[ell]-I-D**kk[0])).rows())
assert Matrix(GF(2),ops).row_space()==W.row_space()and W.nrows()==82
arc=lit((P/'quadratic_cusp_arc_data.txt').read_text());width,pol,prec,pivot,inds,sidx,co,xco,yco=arc;assert width==3 and prec==800;R=PolynomialRing(GF(2),'a');F=GF(4,'a',modulus=R(pol));a=F.gen();PS=PowerSeriesRing(F,'q',default_prec=5);q=PS.gen();z=[PS([sum(F(c)*a**i for i,c in enumerate(v))for v in row[:5]]).add_bigoh(5)for row in co]
assert z[pivot-1]==1 and any(v[1]for v in z)
def norm(zz):
 k=next(i for i,v in enumerate(zz)if v[0]);return[v/zz[k]for v in zz]
def act(M,z):return[sum(int(M[i,j])*z[j]for j in range(85))for i in range(85)]
assert norm(act(D**9,z))==norm([PS([v[n]**2 for n in range(5)]).add_bigoh(5)for v in z])
def bits(v):return sum(int(c)<<i for i,c in enumerate(v))
rats=[];quads=[];qs=set()
for k in range(18):
 B=W*(D**k)*Q[:,1:10];rats.append([bits(B.column(j))for j in range(9 if k==0 else 4)])
for k in range(9):
 zk=act(D**k,z);normz=norm(zk);ky=min(tuple(int(c)for v in normz for c in v[0]._vector_()),tuple(int(c)for v in normz for c in(v[0]**2)._vector_()));assert ky not in qs;qs.add(ky)
 vals=act(W,zk)if False else[sum(int(W[i,j])*zk[j]for j in range(85))for i in range(W.nrows())];cols=[]
 for n in range(4):
  vs=[v[n]for v in vals];cols.extend([bits([int(v.trace())for v in vs]),bits([int((a*v).trace())for v in vs])])
 quads.append(cols)
# Complete prefix enumeration. Every odd-degree cusp divisor has a rational
# support point. Move one with maximal rational multiplicity to cusp0.
def insert(b,v):
 while v:
  i=v.bit_length()-1
  if i in b:v^=b[i]
  else:b[i]=v;return True
 return False
class Bad(Exception):pass
records=[];total=0
for n0 in range(1,10):
 base={}
 for v in rats[0][:n0]:assert insert(base,v)
 groups=rats[1:]+quads;costs=[1]*17+[2]*9;limits=[min(n0,4)]*17+[4]*9;count=[0];bad=[]
 def visit(j,left,b,chosen):
  if j==26:
   if left==0:
    if len(b)!=9:raise Bad((n0,chosen,len(b)))
    count[0]+=1
   return
  if left>sum(costs[k]*limits[k]for k in range(j,26))or(j>=17 and left%2):return
  cost=costs[j];lim=min(limits[j],left//cost);cap=sum(costs[k]*limits[k]for k in range(j+1,26));cur=b.copy()
  for m in range(lim+1):
   if left-m*cost<=cap:
    visit(j+1,left-m*cost,cur,chosen+[m])
   if m<lim:
    ok=True
    for v in groups[j][m*cost:(m+1)*cost]:ok=insert(cur,v)and ok
    if not ok:pass
 try:visit(0,9-n0,base,[])
 except Bad as exc:
  fail=dict(n0=n0,full_choices=exc.args[0][1],rank=exc.args[0][2]);(P/'all_cusp_FI_failure.json').write_text(json.dumps(fail,indent=2)+'\n');print('C57_ALL_CUSP_FI_DEPENDENCY',fail,flush=True);raise
 # Independent generating-function coefficient count.
 dp=[0]*10;dp[0]=1
 for c,lim in zip(costs,limits):
  nd=[0]*10
  for v in range(10):
   for m in range(lim+1):
    if v+c*m<10:nd[v+c*m]+=dp[v]
  dp=nd
 assert count[0]==dp[9-n0];total+=count[0];records.append(dict(max_rational_multiplicity=n0,configurations=count[0]));print('C57_ALL_CUSP_FI_PREFIX_PASS',n0,count[0],'seconds',time.monotonic()-t0,flush=True)
out=dict(complete=True,level=57,prime=2,genus=85,annihilating_differential_dimension=82,rational_cusps=18,quadratic_cusps=9,quadratic_tate_arc_verified=True,all_normalized_configurations=total,by_maximum_rational_multiplicity=records,all_FI_ranks_equal9=True)
if '--readonly'in sys.argv:assert json.loads((P/'all_cusp_FI_sage_certificate.json').read_text())==out
else:(P/'all_cusp_FI_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C57_ALL_GLOBAL_CUSP_DIVISORS_FI_INDEPENDENT_SAGE_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
