#!/usr/bin/env python3
"""Exact independent torsion-family/unit-jet and finite-rank replay."""
from pathlib import Path
import json,time,hashlib,sys,os
os.environ["OPENBLAS_NUM_THREADS"]="1";os.environ["OMP_NUM_THREADS"]="1"
from sage.all import GF,PolynomialRing,matrix,EllipticCurve
from unit_arithmetic93 import *
H=Path(__file__).resolve().parents[1]/'data/N93_source_units';t=time.monotonic();a=(H/'canonical_unit_data.m').read_text();b=(H/'elliptic_hecke5_v2_data.m').read_text();RX=matbits(a,'UnitXRelations');RY=matbits(a,'UnitYRelations');ss=lit(b,'SourceRecords');tt=lit(b,'SourceHeckeRecords');assert len(ss)==len(tt)==1;cache={};count=0
for sid,(s,ts)in enumerate(zip(ss,tt),1):
 assert len(ts)==3 and sorted(d for e,d,m,r in ts)==[18,18,18]and all(m==1 for e,d,m,r in ts)
 for rec in[s]+[r for e,d,m,r in ts]:
  v=verify_record(rec,RX,RY);F=v['F'];p=F.p
  if p not in cache:
   R=PolynomialRing(GF(2),'a');mp=R(v['field_poly']);assert mp.is_irreducible();cache[p]=GF(2**F.d,'a',modulus=mp)
  FF=cache[p];ec=lambda z:FF(R([(z>>j)&1 for j in range(F.d)]));E=EllipticCurve(FF,list(map(ec,v['ainv'])));P=E(list(map(ec,v['point'])));assert P!=E(0)and 93*P==E(0)and 31*P!=E(0)and 3*P!=E(0);u,w=v['units'];cv={n:FF(R([(n>>j)&1 for j in range(F.d)]))for n in {0,1,u,u^1,w,w^1}};M=matrix(FF,[[cv[((row>>i)&1)^(u if(row>>(261+i))&1 else 0)]for i in range(261)]for row in RX]+[[cv[((row>>i)&1)^(w if(row>>(261+i))&1 else 0)]for i in range(261)]for row in RY]);assert M.rank()==260;count+=1
 print('NONIC93_JET_REPLAY_SOURCE',sid,'all_four_families_valid=true',flush=True)
# Independent six-kernel5-isogeny transport, compared by Frobenius-invariant Tate-unit keys.
from collections import Counter
src=decode_record(ss[0]);tars=[decode_record(rec)for e,d,m,rec in tt[0]]
RS=PolynomialRing(GF(2),'v');basepoly=RS(src['field_poly']);targetpoly=RS(tars[0]['field_poly']);L=GF(2**18,'c',modulus=targetpoly)
root=basepoly.change_ring(L).roots()[0][0];SB=cache[src['F'].p];emb=SB.hom([root],L)
FL=Field(bits(tars[0]['field_poly']))
def bcode(z):return bits([int(c)for c in z.polynomial().list()])
def units_general(aa,x,y):
 a1,a2,a3,a4,a6=aa;a3p=2*y+a3+a1*x;a2p=3*x+a2;a4p=3*x*x+2*x*a2+a4-a1*y;a1p=2*a4p/a3p+a1;a2pp=a2p-(a4p/a3p)**2-a1*a4p/a3p;b=-a2pp**3/a3p**2;c=-(a1p*a2pp-a3p)/a3p;r=b/c;s=c*c/(b-c);t=r*s-2*r+1;return[(s-r)/t,t/(s*s-s-r+1)]
E0=EllipticCurve(SB,[SB(RS([(v>>j)&1 for j in range(9)]))for v in src['ainv']]);P0=E0([SB(RS([(v>>j)&1 for j in range(9)]))for v in src['point']]);EL=EllipticCurve(L,[emb(v)for v in E0.a_invariants()]);PL=EL([emb(P0[0]),emb(P0[1])])
phis=EL.isogenies_prime_degree(5);assert len(phis)==6
observed=Counter()
for phi in phis:
 EP=phi.codomain();QP=phi(PL);assert 93*QP==EP(0)and 31*QP!=EP(0)and 3*QP!=EP(0)
 uv=units_general(EP.a_invariants(),QP[0],QP[1]);key=FL.key([1]+[bcode(v)for v in uv],18);observed[key]+=1
expected=Counter()
for (e,d,m,rec),v in zip(tt[0],tars):assert e==2 and d==18 and m==1;expected[v['F'].key([1]+v['units'],18)]+=e
assert observed==expected and sum(observed.values())==6
out={'status':'PASS','sources':1,'ordinary93_torsion_families':count,'precision':5,'unit_relations':[162,157],'all_finite_unit_ranks':260,'all_psi93_derivatives_nonzero':True,'all_canonical_tangents_nonzero':True,'independent_Sage_5_isogenies':6,'actual_level_transport_Frobenius_keys_match':True,'native_data_sha256':hashlib.sha256(b.encode()).hexdigest()}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(out,indent=2)+'\n')
print('NONIC93_ELLIPTIC_HECKE_JETS_INDEPENDENT_PASS',out)
