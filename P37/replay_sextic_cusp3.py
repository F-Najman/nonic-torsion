#!/usr/bin/env python3
"""One-basis proof of all1,140 degree6+cuspidal3 restrictions for C37."""
from pathlib import Path
import re,ast,json,sys,hashlib,time
from itertools import combinations_with_replacement
from collections import Counter
H=Path(__file__).resolve().parent/'data';P=H/'sextic_cusp3_jets.m';s=P.read_text();t0=time.monotonic()
assert s.count('NONIC37_SEXTIC_CUSP3_JET_DATA_READY')==1
assert re.search(r'N:=37;D:=9;gx:=40;udim:=38;',s)
assert ast.literal_eval(re.search(r'ds:=(.*?);',s)[1])==[6]+[1]*18
assert ast.literal_eval(re.search(r'caps:=(.*?);',s)[1])==[1]+[3]*18

def parse(pat):
 out={}
 for m in re.finditer(pat+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S):
  v=m.groups();key=tuple(map(int,v[:-3]));r,c=map(int,v[-3:-1]);a=ast.literal_eval(v[-1]);assert len(a)==r*c and set(a)<={0,1}and key not in out
  out[key]=(r,c,[sum(a[i*c+j]<<j for j in range(c))for i in range(r)],[sum(a[i*c+j]<<i for i in range(r))for j in range(c)])
 return out

def basis(values):
 b={}
 for v in values:
  while v:
   k=v.bit_length()-1
   if k in b:v^=b[k]
   else:b[k]=v;break
 return b

def rem(v,b):
 while v:
  k=v.bit_length()-1
  if k not in b:return v
  v^=b[k]
 return 0

def act(v,M):
 a=0
 while v:
  k=v&-v;a^=M[k.bit_length()-1];v^=k
 return a
CM=parse('CartierMatrix')[()];MD=parse('DiamondMatrix')[()];U=parse('ImageBasis')[()]
assert CM[:2]==MD[:2]==(40,40)and U[:2]==(38,40)
I=[1<<i for i in range(40)];A=[MD[2][i]^I[i]for i in range(40)];ab=basis(A);ub=basis(U[2]);assert len(ab)==len(ub)==38 and all(rem(v,ub)==0 for v in ab.values())
assert all(act(MD[2][i],CM[2])==act(CM[2][i],MD[2])for i in range(40))
power=I
for _ in range(18):power=[act(v,MD[2])for v in power]
assert power==I
BJ=parse(r'basisjets\[(\d+)\]');LJ=parse(r'localjets\[<(\d+),(\d+)>\]');assert len(BJ)==19 and len(LJ)==55
for (i,),block in BJ.items():
 width=6 if i==1 else 3;assert block[:2]==(40,width)
 expected=[sum(((v&col).bit_count()&1)<<k for k,v in enumerate(U[2]))for col in block[3]]
 for m in range(1,2 if i==1 else 4):
  z=LJ[i,m];c=6 if i==1 else m;assert z[:2]==(38,c)and z[3]==expected[:c]
count=0;bad=[]
for v in combinations_with_replacement(range(2,20),3):
 cs=LJ[1,1][3].copy()
 for i,m in Counter(v).items():cs+=LJ[i,m][3]
 rank=len(basis(cs));assert len(cs)==9
 if rank<9:bad.append({'rational_cusp_multiset':v,'rank':rank})
 count+=1
assert count==1140
res={'level':37,'degree':9,'source_genus':40,'alpha':'delta2-1','cotangent_image_dimension':38,'all_rational_cusps':18,'one_NC6_diamond_representative':True,'tested_divisors':count,'exceptions_count':len(bad),'exceptions':bad,'geometry_sha256':hashlib.sha256(P.read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t0,'status':'PASS'if not bad else'FAIL'}
if '--readonly'in sys.argv:
 saved=json.loads((H/'sextic_cusp3_verified.json').read_text());assert all(saved[k]==v for k,v in res.items()if k!='elapsed_seconds')
else:(H/'sextic_cusp3_verified.json').write_text(json.dumps(res,indent=2)+'\n')
print('NONIC37_SEXTIC_CUSP3_REPLAY_COMPLETE',json.dumps(res,sort_keys=True));assert not bad
