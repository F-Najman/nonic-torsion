#!/usr/bin/env python3
"""All C65 cyclic3-kernel isogenies, using kernel polynomials over finite fields."""
from pathlib import Path
import json,time,sys,collections,math,os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from sage.all import GF,PolynomialRing,EllipticCurve,divisors
from unit_arithmetic import *
H=Path(__file__).resolve().parents[1]/'data'/'composite65';t=time.monotonic();b=(H/'elliptic_hecke_jets_stage2_data.m').read_text();ss=lit(b,'SourceRecords');tt=lit(b,'SourceHeckeRecords');fds=lit(b,'SourceDegrees');results=[]
def units(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P[0],P[1];aa3=a3+a1*x;aa2=x+a2;aa4=x*x+a4+a1*y;w=aa4/aa3;aaa2=aa2+w*w+a1*w;b=aaa2**3/aa3**2;c=(a1*aaa2+aa3)/aa3;r=b/c;s=c*c/(b+c);v=r*s+1;return((s+r)/v,v/(s*s+s+r+1))
for sid,(f,s,trs)in enumerate(zip(fds,ss,tt),1):
 v=decode_record(s);dd=math.lcm(*[d for e,d,m,r in trs]);L=GF(2**dd,'z');R=PolynomialRing(L,'X');X=R.gen();FB=Field(bits([int(c)for c in L.modulus()]));enc=lambda a:bits([int(c)for c in a.polynomial()]);root=R(v['field_poly']).roots(multiplicities=False)[0]
 def embed(a):return sum(L((a>>j)&1)*root**j for j in range(f))
 E=EllipticCurve(L,list(map(embed,v['ainv'])));P=E(list(map(embed,v['point'])));assert 65*P==E(0)and 5*P!=E(0)and 13*P!=E(0);ks=(X**4+X**3+E.a6()).roots(multiplicities=False);assert len(ks)==4;actual=collections.Counter()
 for x in ks:
  # A cyclic3-kernel is determined by its monic linear x-polynomial;
  # its individual nonzero points need not be rational over L.
  phi=E.isogeny(X-x);assert phi.degree()==3;ET=phi.codomain();PT=phi(P);assert 65*PT==ET(0)and 5*PT!=ET(0)and 13*PT!=ET(0);xy=units(ET,PT);d=min(d for d in divisors(dd)if all(FB.frob(enc(a),int(d))==enc(a)for a in xy));key=FB.key([1]+[enc(a)for a in xy],int(d));actual[key]+=1
 expected=collections.Counter()
 for e,d,m,rec in trs:
  w=decode_record(rec);expected[w['F'].key([1]+w['units'],d)]+=e
 assert actual==expected;results.append({'source':sid,'source_degree':f,'kernel_polynomials':4,'target_degrees':[d for e,d,m,r in trs]});print('NONIC65_SAGE_KERNEL_ISOGENY_SOURCE_PASS',sid,flush=True)
out={'status':'PASS','sources':results,'independent_actual3_isogenies':28,'all65_level_transports_exact':True};p=H/'hecke_transport_independent_verified.json'
if '--readonly'in sys.argv:assert json.loads(p.read_text())==out
else:p.write_text(json.dumps(out,indent=2)+'\n')
print('NONIC65_SAGE_KERNEL_ISOGENY_TRANSPORT_PASS elapsed',time.monotonic()-t)
