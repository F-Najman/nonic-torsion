from pathlib import Path
import ast,re,json,time,sys
from sage.all import GF,Matrix,PolynomialRing,PowerSeriesRing
P=Path(__file__).resolve().parents[1]/'data'/'allbad';t0=time.monotonic();g=121

def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']'))
def mat(s,name):
 m=re.search(name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,name
 return Matrix(GF(2),int(m[1]),int(m[2]),lit(m[3]))
s=(P.parent/'source_units'/'canonical_unit_data.m').read_text();Q=mat(s,'qcoeffs');D1=mat(s,'Diamond1Matrix');D2=mat(s,'Diamond2Matrix');GG=[D1**a*D2**b for a in range(2)for b in range(12)];assert len({tuple(M.list())for M in GG})==24
s=(P/'annihilator_data.m').read_text();W=mat(s,'AnnihilatingDifferentials');assert W.nrows()==106 and(W*D1).row_space()==W.row_space()and(W*D2).row_space()==W.row_space()
ns=lit(re.search(r'HeckeIndices:=(\[.*?\]);',s,re.S)[1]);aq=lit(re.search(r'PositiveFactorEigenvalues:=(\[.*?\]);',s,re.S)[1]);raw=lit(re.search(r'HeckeMatrices:=(\[.*?\]);AnnihilatorMetadata',s,re.S)[1]);T={n:Matrix(GF(2),g,g,r)for n,r in zip(ns,raw)};I=Matrix.identity(GF(2),g);ops=[]
for n,a in zip(ns,aq):
 for ell in[3,7,11,17]:
  kk=[(u,v)for u in range(2)for v in range(12)if(2**u*3**v)%65 in[ell%65,(-ell)%65]];assert len(kk)==1;u,v=kk[0];ops+=list(((T[n]-GF(2)(a)*I)*(T[ell]-I-D1**u*D2**v)).rows())
assert Matrix(GF(2),ops).row_space()==W.row_space()
width,pol,prec,pivot,inds,sidx,co,xco,yco=lit((P/'quartic_cusp_arc_data.txt').read_text());assert width==5 and prec==800;F=GF(16,'a',modulus=PolynomialRing(GF(2),'x')(pol));a=F.gen();PS=PowerSeriesRing(F,'q',default_prec=5);q=PS.gen();z=[PS([sum(F(c)*a**i for i,c in enumerate(v))for v in row[:5]]).add_bigoh(5)for row in co];assert(z[pivot-1]-1).valuation()>=5 and any(v[1]for v in z)
def norm(zz):
 k=next(i for i,v in enumerate(zz)if v[0]);return[v/zz[k]for v in zz]
def act(M,z):return[sum(int(M[i,j])*z[j]for j in range(g))for i in range(M.nrows())]
def key(zz):
 zz=norm(zz);return min(tuple(int(c)for v in zz for c in(v[0]**(2**k))._vector_())for k in range(4))
seen=set();quartic_arcs=[]
for M in GG:
 zz=act(M,z);ky=key(zz)
 if ky not in seen:seen.add(ky);quartic_arcs.append(zz)
assert len(quartic_arcs)==6
# All four coefficient conjugates of the base arc occur in the full diamond orbit.
conjugates=[norm([PS([v[n]**(2**k)for n in range(5)]).add_bigoh(5)for v in z])for k in range(4)];got=[norm(act(M,z))for M in GG if key(act(M,z))==key(z)];assert len(got)==4 and all(any(all((u-v).valuation()>=5 for u,v in zip(zz,cc))for zz in got)for cc in conjugates)
def bits(v):return sum(int(c)<<i for i,c in enumerate(v))
rats=[];quarts=[]
for k,M in enumerate(GG):
 B=W*M*Q[:,1:10];rats.append([bits(B.column(j))for j in range(9 if k==0 else 4)])
for zz in quartic_arcs:
 vals=act(W,zz);cols=[]
 for n in range(2):
  vs=[v[n]for v in vals]
  for j in range(4):cols.append(bits([int((a**j*v).trace())for v in vs]))
 quarts.append(cols)
def insert(b,v):
 while v:
  k=v.bit_length()-1
  if k in b:v^=b[k]
  else:b[k]=v;return True
 return False
class Bad(Exception):pass
records=[];total=0
for n0 in range(1,10):
 base={}
 for v in rats[0][:n0]:assert insert(base,v),('base',n0)
 groups=rats[1:]+quarts;costs=[1]*23+[4]*6;limits=[min(n0,4)]*23+[2]*6;count=[0]
 def visit(j,left,b,chosen):
  if j==29:
   if left==0:
    if len(b)!=9:raise Bad((n0,chosen,len(b)))
    count[0]+=1
   return
  if left>sum(costs[k]*limits[k]for k in range(j,29))or(j>=23 and left%4):return
  c=costs[j];lim=min(limits[j],left//c);cap=sum(costs[k]*limits[k]for k in range(j+1,29));cur=b.copy()
  for m in range(lim+1):
   if left-m*c<=cap:visit(j+1,left-m*c,cur,chosen+[m])
   if m<lim:
    for v in groups[j][m*c:(m+1)*c]:insert(cur,v)
 try:visit(0,9-n0,base,[])
 except Bad as exc:
  fail=dict(max_rational_multiplicity=n0,choices=exc.args[0][1],rank=exc.args[0][2]);(P/'all_cusp_FI_failure.json').write_text(json.dumps(fail,indent=2)+'\n');print('C65_ALL_CUSP_FI_DEPENDENCY',fail,flush=True);raise
 dp=[0]*10;dp[0]=1
 for c,lim in zip(costs,limits):
  nd=[0]*10
  for v in range(10):
   for m in range(lim+1):
    if v+c*m<10:nd[v+c*m]+=dp[v]
  dp=nd
 assert count[0]==dp[9-n0];total+=count[0];records.append(dict(max_rational_multiplicity=n0,configurations=count[0]));print('C65_ALL_CUSP_FI_PREFIX_PASS',n0,count[0],round(time.monotonic()-t0,2),flush=True)
out=dict(complete=True,level=65,prime=2,genus=121,annihilating_differential_dimension=106,rational_cusps=24,quartic_cusps=6,genuine_quartic_Tate_arc=True,all_normalized_configurations=total,by_maximum_rational_multiplicity=records,all_FI_ranks_equal9=True)
if '--readonly'in sys.argv:assert json.loads((P/'all_cusp_FI_sage_certificate.json').read_text())==out
else:(P/'all_cusp_FI_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C65_ALL_GLOBAL_CUSP_DIVISORS_FI_INDEPENDENT_SAGE_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
