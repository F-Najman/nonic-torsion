#!/usr/bin/env python3
"""Exact norm relation by finite-place CRT and rational reconstruction.
No Magma is launched. The geometric bidegree bound is (31,480), so
972 distinct good base-prime degrees determine each rational coefficient
of height<=480 uniquely (972>2*480).
"""
from pathlib import Path
import ast,re,time,json,sys,hashlib
from sage.all import GF,PolynomialRing,matrix
H=Path(__file__).resolve().parent/'data';t0=time.monotonic();REPLAY='--readonly'in sys.argv;s=(H/'finite_j_model_snapshot.m').read_text();BOUND=480;N=31

def lit(name):return ast.literal_eval(re.search(re.escape(name)+r':=(.*?);',s,re.S)[1].replace('<','(').replace('>',')'))
def bits(p):return sum(int(a)<<i for i,a in enumerate(p))
def pdiv(a,b):
 assert b;q=0;db=b.bit_length()-1
 while a and a.bit_length()-1>=db:
  k=a.bit_length()-1-db;q^=1<<k;a^=b<<k
 return q,a

def rem(a,b):return pdiv(a,b)[1]
def mul(a,b):
 if a.bit_count()<b.bit_count():a,b=b,a
 z=0
 while b:
  k=b&-b;z^=a<<(k.bit_length()-1);b^=k
 return z

def gcd(a,b):
 while b:a,b=b,rem(a,b)
 return a

def exact(a,b):q,r=pdiv(a,b);assert r==0;return q

def inv(a,m):
 r0,r1=m,rem(a,m);s0,s1=0,1
 while r1:q,r=pdiv(r0,r1);r0,r1=r1,r;s0,s1=s1,s0^mul(q,s1)
 assert r0==1;return rem(s0,m)

def reconstruct(r,m):
 if r==0:return 0,1
 r0,r1=m,r;s0,s1=0,1
 while r1.bit_length()-1>BOUND:
  q,r2=pdiv(r0,r1);r0,r1=r1,r2;s0,s1=s1,s0^mul(q,s1)
 assert s1 and s1.bit_length()-1<=BOUND
 gg=gcd(r1,s1);a,b=exact(r1,gg),exact(s1,gg)
 assert gcd(b,m)==1 and rem(a^mul(r,b),m)==0
 assert max(a.bit_length()-1,b.bit_length()-1)<=BOUND
 return a,b

fc=[(bits(a),bits(b))for a,b in lit('AFPolynomialData')];jc=[(bits(a),bits(b))for a,b in lit('JFunctionData')];assert len(fc)==N+1 and len(jc)==N
R=PolynomialRing(GF(2),'t');records=[];skipped=[];M=1;values=[0]*(N+1)
for mi in range((1<<12)+1,1<<13,2):
 mp=R([(mi>>i)&1 for i in range(13)])
 if not mp.is_irreducible():continue
 if any(rem(b,mi)==0 for a,b in fc+jc)or rem(fc[-1][0],mi)==0:skipped.append(mi);continue
 F=GF(2**12,'a',modulus=mp)
 def elem(a):
  a=rem(a,mi);return F(R([(a>>i)&1 for i in range(12)]))
 def rational(c):return elem(c[0])/elem(c[1])
 f=[rational(c)for c in fc];j=[rational(c)for c in jc];lead=f[-1];assert lead;f=[c/lead for c in f];cols=[];v=j
 for _ in range(N):
  cols.append(v);a=v[-1];v=[F(0)]+v[:-1]
  if a:v=[x+a*y for x,y in zip(v,f[:-1])]
 A=matrix(F,cols).transpose();cp=A.charpoly();assert cp.degree()==N and cp[N]==1 and cp[N-1]==A.trace()
 residues=[bits(c.polynomial().list())for c in cp.list()];assert len(residues)==N+1
 imin=inv(rem(M,mi),mi)
 for i,a in enumerate(residues):delta=rem(mul(a^rem(values[i],mi),imin),mi);values[i]^=mul(M,delta)
 M=mul(M,mi);records.append({'modulus':mi,'characteristic_polynomial':residues})
 if len(records)%10==0:print('NONIC87_MODULAR_NORM_PROGRESS primes',len(records),'degree',M.bit_length()-1,'elapsed',time.monotonic()-t0,flush=True)
 if M.bit_length()-1>2*BOUND:break
assert len(records)==81 and M.bit_length()-1==972
fractions=[reconstruct(a,M)for a in values];assert fractions[-1]==(1,1)
D=1
for a,b in fractions:D=mul(exact(D,gcd(D,b)),b)
cc=[mul(a,exact(D,b))for a,b in fractions];content=0
for a in cc:content=gcd(content,a)
cc=[exact(a,content)for a in cc];assert len(cc)==32 and max(a.bit_length()-1 for a in cc)==480 and cc[-1]
# Check every reconstructed coefficient against every exact residue again.
for rec in records:
 m=rec['modulus']
 for (a,b),v in zip(fractions,rec['characteristic_polynomial']):assert rem(a^mul(v,b),m)==0 and gcd(b,m)==1
out={'status':'PASS','relation':'field norm characteristic polynomial; may be a power of the minimal polynomial','base_degree':31,'j_degree':480,'rational_height_bound':480,'base_primes':records,'skipped_bad_primes':skipped,'CRT_modulus':M,'CRT_modulus_degree':972,'rational_coefficients':[[a,b]for a,b in fractions],'primitive_relation_coefficients':cc,'plane_model_sha256':hashlib.sha256(s.encode()).hexdigest(),'elapsed_seconds':time.monotonic()-t0}
f=H/'modular_jbase_norm_verified.json'
if REPLAY:
 old=json.loads(f.read_text())
 for k in out:
  if k!='elapsed_seconds':assert old[k]==out[k],k
else:
 f.write_text(json.dumps(out,indent=2)+'\n')
 cofs=[[(a>>i)&1 for i in range(max(1,a.bit_length()))]for a in cc]
 (H/'jbase_relation_modular_data.m').write_text('JBaseDegreeJ:=31;JBaseDegreeParameter:=480;JBasePolynomialCoefficients:='+str(cofs)+';\nprint "NONIC87_MODULAR_JBASE_RELATION_READY";\n')
print('NONIC87_MODULAR_NORM_CRT_CERTIFICATE_PASS primes81 modulusdegree972 height480 bidegree31x480 elapsed',time.monotonic()-t0,flush=True)
