#!/usr/bin/env python3
"""Exact independent integer replay; no native class-group computation."""
from pathlib import Path
import ast,json
P=Path(__file__).resolve().parent/'data'

def load(name):
 s=(P/name).read_text();D={}
 for n in ast.parse(s).body:
  assert isinstance(n,ast.Assign) and len(n.targets)==1 and isinstance(n.targets[0],ast.Name)
  D[n.targets[0].id]=ast.literal_eval(n.value)
 return D

def residue(v,mods):
 assert len(v)==4 and v[3]==0
 return tuple(v[i]%mods[i] for i in range(3))

def check_context(D):
 assert D['relations2']==[[3,0,0,0],[0,3,0,0],[0,0,52497,0]]
 assert D['relations5']==[[3,0,0,0],[0,3,0,0],[0,0,342175446,0]]
 # Exact HNF computation independently proves the membership shortcut.
 from sympy import Matrix
 from sympy.matrices.normalforms import hermite_normal_form
 H=hermite_normal_form(Matrix([r[:3] for r in D['relations5']+D['cusp5']]).T)
 assert H==Matrix.diag(1,1,6518)
 return {residue(v,[3,3,52497]) for v in D['candidate2']}

def replay_rows(D,ordinary=None):
 candidate=check_context(D);count=members=passes=0;survivors=[]
 for label,aa,mem,pas,im,pre in D['tests']:
  assert len(aa)==4 and aa[3]==0
  m=(aa[2]%6518==0);assert int(m)==mem
  pp=False
  if m:
   members+=1;assert len(pre)==14
   a=[sum(v*c[j] for v,c in zip(pre,D['cusp5'])) for j in range(4)]
   b=[sum(v*c[j] for v,c in zip(pre,D['cusp2'])) for j in range(4)]
   assert residue(a,[3,3,342175446])==residue(aa,[3,3,342175446])
   assert residue(b,[3,3,52497])==residue(im,[3,3,52497])
   pp=residue(im,[3,3,52497]) in candidate
  else:assert im==pre==[]
  assert int(pp)==pas;count+=1;passes+=pp
  if pp and (ordinary is None or residue(im,[3,3,52497]) in ordinary):survivors.append(label)
 return {'tests':count,'members':members,'passes':passes,'survivors':survivors}

if __name__=='__main__':
 A=load('filtered9_classes_shard0.py');assert A['complete'];ord2=check_context(A);assert len(ord2)==2412
 L=load('lower_testclasses.py');assert L['complete'] and L['cusp2']==A['cusp2'];lower=replay_rows(L,ord2)
 assert lower['tests']==278876 and lower['members']==641 and lower['passes']==10
 assert lower['survivors']==[[3,6,17,1531],[4,5,102,364],[4,5,23,190]]
 octic=[];old_indices=set()
 for sh in range(4):
  D=load(f'degree8/old_good8_classes_shard{sh}.py');assert D['complete']
  R=replay_rows(D);assert R['members']==0

  for row in D['tests']:
   idx=row[0];assert idx not in old_indices;old_indices.add(idx);assert (idx-1)%4==sh
  octic.append(R)
 assert old_indices==set(range(1,3580))
 import sys
 if '--lower-only' in sys.argv:
  report={'lower':lower,'octic_total':sum(r['tests'] for r in octic),'octic_H5_members':0,'complete':True}
  print('NONIC54_LOWER_OCTIC_INTEGER_REPLAY '+json.dumps(report));raise SystemExit
 seen=set();skip=set();pure=[]
 for sh in range(4):
  D=load(f'filtered9_classes_shard{sh}.py');assert D['complete'] and D['cusp2']==A['cusp2'];assert check_context(D)==ord2
  R=replay_rows(D);pure.append(R)
  for row in D['tests']:
   idx=row[0];assert idx not in seen and (idx-1)%4==sh;seen.add(idx)
  ss=set(D['quadratic_skip_indices']);assert not ss&skip;assert all((i-1)%4==sh for i in ss);skip|=ss
 Q=json.loads((P/'quadratic_prepass_results.json').read_text())
 assert len(seen)==8048 and len(skip)==8087 and not seen&skip and seen|skip==set(range(1,16136))
 assert skip==set(Q['nonsquare_indices'])
 report={'lower':lower,'octic_total':3579,'octic_H5_members':0,'pure9_queries':8048,'quadratic_exclusions':8087,'pure9_survivors':sorted(i for r in pure for i in r['survivors']),'all_first_projection_classes_replayed':True,'complete':True}
 if '--write-json' in sys.argv:(P/'complete_classes_integer_replay.json').write_text(json.dumps(report,indent=2)+'\n')
 print('NONIC54_COMPLETE_INTEGER_REPLAY '+json.dumps(report))
