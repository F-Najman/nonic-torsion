#!/usr/bin/env python3
"""Complete exact R1 multiplication sieve for actual Hecke target diamond orbits."""
from pathlib import Path
import json,gzip,time,random,sys
from unit_arithmetic import *
H=Path(__file__).resolve().parents[1]/'data/canonical_units';t0=time.monotonic();REPLAY='--readonly'in sys.argv;rng=random.Random(20261026)
a=(H/'canonical_unit_stage1_data.m').read_text();b=(H/'elliptic_hecke_jets_stage2_data.m').read_text();D=matbits(a,'Diamond2Matrix');sr=lit(b,'SourceRecords');tr=lit(b,'SourceHeckeRecords');assert len(sr)==len(tr)==5;T=json.loads((H/'canonical_ring_R1_tensor.json').read_text());assert T['dimensions']==[247,738];TT=[[int(v,16)for v in row]for row in T['product_1_1_into_2']]
class Span:
 def __init__(self,rows=()):
  self.p={}
  for v in rows:self.add(v)
 def add(self,v):
  while v:
   k=(v&-v).bit_length()-1
   if k in self.p:v^=self.p[k]
   else:self.p[k]=v;return True
  return False
 def quotient(self,v):
  for k in sorted(self.p):
   if v>>k&1:v^=self.p[k]
  return v
 def __len__(self):return len(self.p)
def kernel(rows):
 p={};out=[]
 for i,v in enumerate(rows):
  c=1<<i
  while v:
   k=(v&-v).bit_length()-1
   if k in p:a,b=p[k];v^=a;c^=b
   else:p[k]=(v,c);break
  else:out.append(c)
 return out
families={};keys={};SOURCE_ORBITS=[];TARGET_ORBITS=[]
def orbit(rec,degree):
 v=decode_record(rec);F=v['F'];z=v['jets'];assert F.normalize(z)==z;original=z;oo=[]
 for k in range(41):
  key=F.key(z,degree)
  if key not in keys:
   fid=len(keys)+1;keys[key]=fid;families[fid]={'F':F,'jets':z,'degree':degree,'evals':{}}
  fid=keys[key];assert families[fid]['degree']==degree;oo.append(fid);z=transform(F,z,D)
 assert z==original and len(set(oo))==41;return oo
for sid,(s,targets)in enumerate(zip(sr,tr),1):
 SOURCE_ORBITS.append(orbit(s,9));TARGET_ORBITS.append([(mult,orbit(r,deg))for e,deg,mult,r in targets]);assert sum(m*families[o[0]]['degree']for m,o in TARGET_ORBITS[-1])==36
 print('NONIC83_DIAMOND_JET_ORBITS_READY source',sid,'families_so_far',len(families),'elapsed',time.monotonic()-t0,flush=True)
assert len(set(sum(SOURCE_ORBITS,[])))==205

def evaluation(fid,m):
 v=families[fid];assert 1<=m<=5
 if m not in v['evals']:
  F=v['F'];v['evals'][m]=[sum(a[n]<<(F.d*n)for n in range(m))for a in v['jets']]
 return v['evals'][m],v['F'].d*m

def conditions(dd,positive):
 rows=[0]*247;width=0;degree=0;support=[]
 for fid,c in sorted(dd.items()):
  c=c if positive else-c
  if c<=0:continue
  ev,w=evaluation(fid,c);rows=[x|(y<<width)for x,y in zip(rows,ev)];width+=w;degree+=families[fid]['degree']*c;support.append([fid,c])
 return rows,width,degree,support
saved=json.loads(gzip.open(H/'canonical_R1_matrices.json.gz','rt').read())if REPLAY else None;cases=[];matrices=[];unresolved=[];ix=0
for sid,(so,to)in enumerate(zip(SOURCE_ORBITS,TARGET_ORBITS),1):
 for shift in range(41):
  for conv in[1,2]:
   dd={}
   def add(fid,n):dd[fid]=dd.get(fid,0)+n
   for mult,o in to:add(o[(shift+1)%41],mult);add(o[shift],-mult)
   # Canonical D is point diamond2^{-1}; point diamond3 = D^10.
   if conv==1:add(so[11],-1);add(so[10],1);add(so[1],-3);add(so[0],3)
   else:add(so[1],-1);add(so[0],1);add(so[11],-3);add(so[10],3)
   bp,pw,pd,ps=conditions(dd,True);bm,nw,nd,ns=conditions(dd,False);assert pd==nd and pd<=72
   kk=kernel(bp);assert len(kk)>=247-pd;old=saved[ix]if REPLAY else None
   fv=int(old['F'],16)if REPLAY else kk[0];assert fv and xorrows(fv,bp)==0
   W=Span([xorrows(fv,[row[j]for row in TT])for j in range(247)]);assert len(W)==247;free=[j for j in range(738)if j not in W.p];assert len(free)==491;rows=bm[:];hh=[]
   candidates=[int(v,16)for v in old['H']]if REPLAY else[xorrows(rng.getrandbits(len(kk)),kk)for _ in range(16)]
   for hv in candidates:
    if not hv or hv==fv:continue
    assert xorrows(hv,bp)==0;hh.append(hex(hv))
    for i in range(247):
     v=W.quotient(xorrows(hv,TT[i]));rows[i]|=sum(((v>>j)&1)<<k for k,j in enumerate(free))<<(nw+491*(len(hh)-1))
    rank=len(Span(rows))
    if rank==247:break
   rec={'source':sid,'target_diamond_translate':shift,'convention':conv,'positive_degree':pd,'positive_conditions_width':pw,'negative_conditions_width':nw,'positive_support':ps,'negative_support':ns,'rank':rank,'F':hex(fv),'H':hh,'columns':nw+491*len(hh)}
   if rank!=247:assert not REPLAY;unresolved.append(rec)
   else:
    mr={k:rec[k]for k in['source','target_diamond_translate','convention','F','H','columns']}|{'rows':[hex(v)for v in rows]}
    if REPLAY:assert mr==old
    matrices.append(mr);cases.append(rec)
   ix+=1
   if ix%82==0:print('NONIC83_R1_SIEVE_PROGRESS cases',ix,'excluded',len(cases),'unresolved',len(unresolved),'elapsed',time.monotonic()-t0,flush=True)
assert ix==410
if REPLAY:assert len(saved)==410;print('NONIC83_R1_SIEVE_READONLY_PASS cases410 rank247 elapsed',time.monotonic()-t0)
else:
 out={'status':'PARTIAL'if unresolved else'PASS','cases':len(cases),'unresolved':unresolved,'canonical_genus':247,'domain_dimension':247,'product_dimension':738,'residual_quotient_dimension':491,'source_diamond_orbits':5,'source_closed_points':205,'diamond_order':41,'point_diamond3_canonical_exponent':10,'all_support_multiplicities_aggregated':True,'results':cases,'elapsed_seconds':time.monotonic()-t0};(H/'canonical_R1_sieve_verified.json').write_text(json.dumps(out,indent=2)+'\n')
 with gzip.open(H/'canonical_R1_matrices.json.gz','wt')as f:json.dump(matrices,f,separators=(',',':'))
 print('NONIC83_R1_SIEVE_RESULT',len(cases),'excluded',len(unresolved),'unresolved elapsed',time.monotonic()-t0,flush=True)
