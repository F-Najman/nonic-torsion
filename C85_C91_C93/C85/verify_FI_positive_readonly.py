#!/usr/bin/env python3
"""Reconstruct all rational-cusp rows and exhaust the entire relation kernel."""
from pathlib import Path
import json,time,sys,hashlib,collections
from sage.all import GF,matrix
from unit_arithmetic import lit,matbits,xorrows
H=Path(__file__).resolve().parents[1]/'data/C85';t=time.monotonic();a=(H/'canonical_unit_stage1_data.m').read_text();b=(H/'FI_positive_data.m').read_text();Q=matbits(a,'qcoeffs');D2=matbits(a,'Diamond2Matrix');D3=matbits(a,'Diamond3Matrix');T2=matbits(b,'CanonicalHecke2');T3=matbits(b,'CanonicalHecke3');g=225;II=[1<<i for i in range(g)]
def prod(A,B):return[xorrows(v,B)for v in A]
def power(A,n):
 B=II[:]
 while n:
  if n&1:B=prod(B,A)
  A=prod(A,A);n>>=1
 return B
for p,D,T in[(2,D2,T2),(3,D3,T3)]:
 length=1300//p if p==3 else 650
 mask=(1<<length)-1;DQ=prod(D,Q);rhs=[sum((((row>>(p*n))&1)^(((p&1)*((drow>>(n//p))&1))if n%p==0 else 0))<<n for n in range(length))for row,drow in zip(Q,DQ)]
 assert [v&mask for v in prod(T,Q)]==rhs
 assert matrix(GF(2),[[v>>j&1 for j in range(length)]for v in Q]).rank()==g
assert lit(b,'KillPolynomialCoefficients')==[[4,12,12,6,1],[-1,0,4,4,1]]
assert lit(b,'WindingFactorData')==[[[2,4,1],[2,2,1]],[[-1,2,1],[1,2,1]]]
A3=[t^i^d for t,i,d in zip(T3,II,D3)];mt=prod(power(T3,4),A3)
labels=[];rows=[]
for aa in range(2):
 for bb in range(16):
  u=pow(2,aa,85)*pow(3,bb,85)%85;labels.append(min(u,85-u));dd=prod(power(D2,aa),power(D3,bb));fq=prod(prod(mt,dd),Q)
  rows.extend(sum(((v>>j)&1)<<i for i,v in enumerate(fq))for j in range(1,10 if len(labels)==1 else 5))
assert labels==lit(b,'DiamondLabels')and len(rows)==133
saved=(H/'FI_positive_rows.txt').read_text().splitlines();assert saved[0]=='85 32 133 450';assert rows==[int(s[:225][::-1],2)for s in saved[1:]]
piv={};kk=[]
for i,v in enumerate(rows):
 c=1<<i
 while v:
  k=v.bit_length()-1
  if k in piv:w,d=piv[k];v^=w;c^=d
  else:piv[k]=(v,c);break
 else:kk.append(c)
assert len(piv)==122 and len(kk)==11
hist=collections.Counter()
for m in range(1,1<<11):
 v=xorrows(m,kk);assert not xorrows(v,rows)
 ns=[(v&511).bit_length()]+[((v>>(9+4*j))&15).bit_length()for j in range(31)];hist[sum(ns)]+=1
assert min(hist)==32 and sum(hist.values())==2047
assert matrix(GF(2),[[v>>j&1 for j in range(225)]for v in rows]).rank()==122
# Count normalized configurations only for information: the kernel argument covers all at once.
configs={}
for n0 in range(1,10):
 dp=[1]+[0]*9
 for _ in range(31):dp=[sum(dp[s-k]for k in range(min(s,n0,4)+1))for s in range(10)]
 configs[n0]=dp[9-n0]
out={'status':'PASS','rational_cusps':32,'candidate_rows':133,'columns':225,'rank':122,'relation_kernel_dimension':11,'nonzero_kernel_vectors_exhausted':2047,'minimum_cusp_prefix_cost':32,'kernel_prefix_cost_histogram':dict(sorted(hist.items())),'all_degree9_prefix_configurations_independent':True,'normalized_configurations_by_maximum_at_first_cusp':configs,'normalized_configurations_total':sum(configs.values()),'finite_image_polynomial_T3':[4,12,12,6,1],'q_Hecke_relations_verified':True,'independent_M4RI_rank':122,'input_sha256':{n:hashlib.sha256((H/n).read_bytes()).hexdigest()for n in['canonical_unit_stage1_data.m','FI_positive_data.m','FI_positive_rows.txt']}}
if '--readonly'in sys.argv:assert json.loads((H/'FI_positive_independent_verified.json').read_text())==out
else:(H/'FI_positive_independent_verified.json').write_text(json.dumps(out,indent=2)+'\n')
print('NONIC85_FI_POSITIVE_INDEPENDENT_PASS rank122 kernel11 minimum_prefix_cost32 elapsed',time.monotonic()-t)
