#!/usr/bin/env sage -python
"""Independent GF(3) replay of canonical multiplication obstruction witnesses.
Reconstructs restrictions from saved local jets, all selected products from the
exact native tensor, annihilator identities and full-rank obstruction matrices.
It does not infer a global exclusion without the source/Hecke and model audits.
"""
from pathlib import Path
import argparse,ast,json,re,time,os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from sage.all import GF,PolynomialRing,PowerSeriesRing,matrix,vector,zero_matrix,ZZ
ap=argparse.ArgumentParser();ap.add_argument('certificates',type=Path);ap.add_argument('--jets',type=Path,required=True);ap.add_argument('--output',type=Path);a=ap.parse_args()
B=Path(__file__).resolve().parent;F=GF(3);t0=time.monotonic()
def stream(s):
 dec=json.JSONDecoder();i=0
 while i<len(s):
  while i<len(s) and s[i].isspace():i+=1
  if i==len(s):break
  row,i=dec.raw_decode(s,i);yield row
raw=(B/'data'/'multiplication_22_data.m').read_text();assert raw.rstrip().endswith('assert #Tensor22Blocks eq 105;')
exps=[]
for d,n in [(2,105),(4,245)]:
 v=ast.literal_eval(re.search(r'Std'+str(d)+r'Exponents:=(.*?);',raw,re.S)[1]);assert len(v)==n and all(len(e)==36 and sum(e)==d for e in v);exps.append(v)
blocks=[matrix(F,105,245,ast.literal_eval(m[1])) for m in re.finditer(r'Append\(~Tensor22Blocks,Matrix\(GF\(3\),105,245,(\[.*?\])\)\);',raw,re.S)]
assert len(blocks)==105
jets=list(stream(a.jets.read_text()));assert jets[0]['type']=='header' and jets[-1]['type']=='complete';cases=jets[-1]['cases']
assert jets[-1]['source_distinct_keys']==10*cases
J={k:{s:[] for s in ('T','U')} for k in range(1,cases+1)}
for r in jets:
 if r['type']=='jet':J[r['case']][r['side']].append(r)
for k in J:
 for side in J[k]:assert sum(r['degree']*r['multiplicity'] for r in J[k][side])==54
sparse=[[[ (i,e) for i,e in enumerate(row) if e] for row in ee] for ee in exps[:1]]
PF=PolynomialRing(F,'x');x=PF.gen()
def eval_jets(r):
 d,m=r['degree'],r['multiplicity'];poly=PF(r['field_polynomial']);assert poly.degree()==d and poly.is_monic()
 K=GF(3**d,name='a',modulus=poly);aa=K.gen();pp=[aa**i for i in range(d)]
 def dec(v):
  assert len(v)==d and all(c in (0,1,2) for c in v)
  return sum((K(c)*u for c,u in zip(v,pp)),K.zero())
 PS=PowerSeriesRing(K,'t',default_prec=m)
 zs=[PS([dec(v) for v in row],prec=m) for row in r['canonical_jets']];assert len(zs)==36 and any(z[0] for z in zs)
 powers={(i,0):PS.one() for i in range(36)}
 for i,z in enumerate(zs):
  for e in range(1,5):powers[i,e]=z if e==1 else powers[i,e-1]*z
 out=[]
 for ee in sparse:
  rows=[]
  for mon in ee:
   v=PS.one()
   for i,e in mon:v*=powers[i,e]
   row=[]
   for k in range(m):
    coeff=v[k].polynomial().list();assert len(coeff)<=d;row += [int(c) for c in coeff]+[0]*(d-len(coeff))
   rows.append(row)
  out.append(matrix(F,rows))
 return out
cert=list(stream(a.certificates.read_text()));assert cert[0]['type']=='header' and cert[-1]['type']=='complete'
assert cert[0]['cases']==cert[-1]['cases'];expected=cert[-1]['cases'];assert expected<=cases
rows=[r for r in cert if r['type']=='case'];assert len(rows)==expected and [r['case'] for r in rows]==list(range(1,expected+1))
assert cert[-1]['survivors']==[]
for r in rows:
 k=r['case'];T=[eval_jets(z) for z in J[k]['T']];U=[eval_jets(z) for z in J[k]['U']]
 def join(v):
  M=v[0]
  for A in v[1:]:M=M.augment(A)
  return M
 ET2=join([z[0] for z in T]);EU2=join([z[0] for z in U])
 assert ET2==matrix(F,105,54,r['ET2']) and EU2==matrix(F,105,54,r['EU2'])
 f=vector(F,r['F']);h=vector(F,r['H']);assert len(f)==105 and len(h)==105 and f and h and not(f*ET2) and not(h*ET2)
 MF=sum((f[i]*blocks[i] for i in range(105) if f[i]),zero_matrix(F,105,245));MH=matrix(F,[h*A for A in blocks]);Q=matrix(F,245,140,r['quotient_projection']);assert not(MF*Q)
 A=EU2.augment(MH*Q);assert A==matrix(F,105,194,r['matrix']);rank=A.rank();assert rank==r['rank']==105
 print('C44_SAGE_MATRIX_CASE case='+str(k)+' rank=105',flush=True)
report={'level':44,'p':3,'cases':expected,'all_local_restrictions_reconstructed':True,'all_product_annihilator_identities':True,'all_obstruction_ranks':105,'tensor_products':11025,'elapsed_seconds':round(time.monotonic()-t0,3),'complete':True,'global_exclusion_requires_separate_source_and_reduction_proof':True}
if a.output:a.output.write_text(json.dumps(report,indent=2)+'\n')
print('C44_SAGE_MATRIX_REPLAY_PASS '+json.dumps(report),flush=True)
