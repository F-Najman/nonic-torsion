#!/usr/bin/env python3
"""Independent integer replay of nonic54 finite class tests."""
from pathlib import Path
import ast,json,sys
from sympy import Matrix,diag
from sympy.matrices.normalforms import hermite_normal_form
P=Path(__file__).resolve().parent/'data'

def load(path):
 if path.suffix=='.json':return json.loads(path.read_text())
 out={}
 for a in ast.parse(path.read_text()).body:
  assert isinstance(a,ast.Assign) and len(a.targets)==1;name=a.targets[0].id;assert name not in out;out[name]=ast.literal_eval(a.value)
 return out

def context(D):
 def finite(v):assert len(v)==4 and v[3]==0;return v[:3]
 def relation(key):
  rows=D[key];assert len(rows)==3 and all(len(v)==4 and v[3]==0 for v in rows);return Matrix([r[:3] for r in rows])
 R2,R5=relation('relations2'),relation('relations5');C2=[finite(v) for v in D['cusp2']];C5=[finite(v) for v in D['cusp5']];assert len(C2)==len(C5)==14
 H2=hermite_normal_form(Matrix(R2.tolist()+C2).T);H5=hermite_normal_form(Matrix(R5.tolist()+C5).T);HJ=hermite_normal_form(Matrix(diag(R2,R5).tolist()+[a+b for a,b in zip(C2,C5)]).T)
 assert abs(R2.det())==472473 and abs(R5.det())==3079579014
 assert abs(H2.det())==1 and abs(H5.det())==6518
 assert abs(R2.det()*R5.det())//abs(HJ.det())==472473
 def member(A,v):return all(a.q==1 for a in A*Matrix(v))
 R2i,R5i,H5i,HJi=R2.T.inv(),R5.T.inv(),H5.inv(),HJ.inv()
 def residue(v):return tuple(a%1 for a in R2i*Matrix(v))
 candidate={residue(finite(v)) for v in D['candidate2']};assert len(candidate) in (2412,2413)
 return finite,member,residue,candidate,C2,C5,R2i,R5i,H5i,HJi

def replay(D,mode):
 finite,member,residue,candidate,C2,C5,R2i,R5i,H5i,HJi=context(D);seen=set();members=passes=count=0
 for row in D['tests']:
  label,q,m,p,im2,preimage=row;q=finite(q);is_member=member(H5i,q);assert int(is_member)==m and p in(0,1)
  if mode=='pure9':assert label not in seen;seen.add(label)
  passed=False
  if is_member:
   members+=1;im2=finite(im2);assert len(preimage)==14 and member(HJi,im2+q)
   got5=[sum(a*c[j] for a,c in zip(preimage,C5)) for j in range(3)];got2=[sum(a*c[j] for a,c in zip(preimage,C2)) for j in range(3)]
   assert member(R5i,[a-b for a,b in zip(got5,q)]) and member(R2i,[a-b for a,b in zip(got2,im2)])
   passed=residue(im2) in candidate
  else:assert im2==preimage==[]
  assert int(passed)==p;passes+=passed;count+=1
 return {'tests':count,'members':members,'passes':passes,'candidate_classes':len(candidate)},seen

def preflight():
 D=load(P/'preflight_components.json');assert D['complete'];finite,member,residue,candidate,*_=context(D)
 assert len(candidate)==2413 and D['noncuspidal_counts2']==[0,0,0,0,0,9,9,18,91]
 nc=D['noncuspidal_classes2'];cs=D['cusp_complement_classes2'];assert list(map(len,nc))==[9,9,18,91]
 assert [len({residue(finite(v)) for v in row}) for row in cs]==[1,9,48,192]
 got=set()
 for i in range(4):
  for a in nc[i]:
   for b in cs[3-i]:got.add(residue([x+y for x,y in zip(finite(a),finite(b))]))
 assert got==candidate
 print('NONIC54_PREFLIGHT_INTEGER_REPLAY '+json.dumps({'candidate_classes':2413,'all_NC6to9_and_cusp_complements_recomputed':True,'cross_prime_cusp_graph_order':472473,'H5_index':6518,'complete':True}))

if __name__=='__main__':
 if len(sys.argv)==1 or sys.argv[1]=='preflight':preflight()
 else:
  D=load(Path(sys.argv[1]));assert D.get('complete');mode='pure9' if 'shard' in D else 'lower';R,seen=replay(D,mode)
  print('NONIC54_CLASS_INTEGER_REPLAY '+json.dumps(R|{'mode':mode,'source_indices':sorted(seen),'supports_exclusion':R['passes']==0,'complete':True}))
