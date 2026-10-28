#!/usr/bin/env python3
"""Independent binary replay of exact RR restrictions and split Hecke images."""
from pathlib import Path
from itertools import product
import ast,json,re,runpy,sys
P=Path(__file__).resolve().parent
F=runpy.run_path(str(P/'replay_local_fi.py'))
mul,add,piv,contains,apply,equalspans=(F[k] for k in ('mul','add','pivots','contains','apply','equalspans'))
source=Path(sys.argv[1]) if len(sys.argv)>1 else P/'data'/'fi3_data_sharp_direct.m'
txt=source.read_text();assert 'FI_JET_DATA_READY' in txt
G=35

def read(name):
 q=re.search(re.escape(name)+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',txt)
 assert q,name
 n,m=map(int,q.group(1,2));return F['matrix'](ast.literal_eval(q.group(3)),n,m),n,m
C,n,m=read('CartierMatrix');assert n==m==G
D,n,m=read('DiamondMatrix');assert n==m==G
I=[1<<i for i in range(G)];D2=mul(D,D)
assert mul(D2,D)==I and mul(D,C)==mul(C,D)
assert G-len(piv(add(D,I)))==11
T=mul(D,C);T2=mul(T,T);T4=mul(T2,T2)
# beta*(1+delta+delta^2)*F0(T), beta*(2-delta-delta^2)*(T-2delta^2), mod2.
U0=mul(mul(add(T,D),add(I,D,D2)),add(T4,T2,I))
U3=mul(mul(add(T,D),add(D,D2)),T)
U=list(piv(U0+U3).values());assert len(piv(U0))==2 and len(piv(U3))==20 and len(U)==22
old,n,m=read('ImageBasis');assert n in (8,22) and m==G
assert equalspans(old,mul(mul(add(T,D),add(T4,T2,I)),T) if n==8 else U)
recorded_u=n
ds=ast.literal_eval(re.search(r'ds:=(\[.*?\]);',txt).group(1));assert sorted(ds)==[1]*6+[3]*2
raw={};blocks={};rr=0
for i,d in enumerate(ds,1):
 A,n,m=read(f'basisjets[{i}]');assert n==G and m==3
 raw[i]=A
 for a in range(1,3//d+1):
  Q=[r&((1<<(a*d))-1) for r in A]
  W,n,m=read(f'kernelspaces[<{i},{a}>]');assert m==G and len(piv(W))==n
  assert all(apply(w,Q)==0 for w in W) and len(piv(Q))+n==G
  J,n,m=read(f'localjets[<{i},{a}>]');assert n==recorded_u and m==a*d and J==mul(old,Q)
  blocks[i,a]=mul(U,Q);rr+=1
assert rr==20
ranks={};bad=[]
for v in product(*(range(3//d+1) for d in ds)):
 if sum(a*d for a,d in zip(v,ds))!=3:continue
 rows=[0]*22;offset=0
 for i,a in enumerate(v,1):
  if not a:continue
  J=blocks[i,a];rows=[x^(y<<offset) for x,y in zip(rows,J)];offset+=a*ds[i-1]
 assert offset==3;rk=len(piv(rows));ranks[v]=rk
 if rk<3:bad.append({'multiplicities':list(v),'rank':rk})
assert len(ranks)==58
report={'genus':35,'finite_operators':2,'U0_dimension':2,'U3_dimension':20,'combined_U_dimension':22,'exact_RR_blocks':rr,'all_cusp_configurations':58,'exceptions':bad,'supports_exclusion':not bad,'complete':True}
print('NONIC133_SPLIT_BINARY_REPLAY '+json.dumps(report))
if '--write-report' in sys.argv:
 (P/'data'/'split_fi_replay.json').write_text(json.dumps(report|{'ranks':[{'multiplicities':list(v),'rank':r} for v,r in ranks.items()]},indent=2)+'\n')
