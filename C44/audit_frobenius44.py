#!/usr/bin/env python3
"""Replay Eichler--Shimura/Newton counts from rational36x36 Hecke matrices."""
from pathlib import Path
import json,sys
from sympy import Matrix,Rational,eye,zeros
p=Path(__file__).resolve().parent;d=json.loads((p/'data'/'N44_frobenius3_source.json').read_text());g=36
T=Matrix(g,g,list(map(lambda v:Rational(*v),d['T3_entries'])));D=Matrix(g,g,list(map(lambda v:Rational(*v),d['Diamond3_entries'])));assert T*D==D*T
J=int((eye(g)-T+3*D).det());assert J==d['local_J_order']==10863413200701600000
S0=2*eye(g);S1=T;tr=[];closed=[]
for n in range(1,10):
 S=S1 if n==1 else T*S1-3*D*S0
 if n>1:S0,S1=S1,S
 t=S.trace();assert t.q==1;tr.append(int(t));points=3**n+1-int(t);c=points-sum(e*closed[e-1] for e in range(1,n) if n%e==0);assert c%n==0;closed.append(c//n)
nc=closed.copy()
for n,c in [(1,15),(2,5),(5,3)]:nc[n-1]-=c
assert closed==d['closed_place_counts']==[15,5,0,0,63,160,240,920,2260]
assert nc==d['noncuspidal_closed_counts']==[0,0,0,0,60,160,240,920,2260]
r={'level':44,'p':3,'local_J_order':J,'Frobenius_traces':tr,'closed_place_counts':closed,'noncuspidal_closed_counts':nc,'independent_rational_matrix_replay':True,'complete':True}
if '--write' in sys.argv:(p/'data'/'frobenius_audit44.json').write_text(json.dumps(r,indent=2)+'\n')
print('C44_FROBENIUS_AUDIT_PASS '+json.dumps(r))
