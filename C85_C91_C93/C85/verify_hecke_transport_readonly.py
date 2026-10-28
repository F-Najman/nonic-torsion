#!/usr/bin/env python3
"""Independent Sage Velu construction of every degree-three target with level."""
from pathlib import Path
import json,time,sys,collections
from sage.all import GF,PolynomialRing,EllipticCurve
from unit_arithmetic import *
H=Path(__file__).resolve().parents[1]/'data/C85';t=time.monotonic();b=(H/'elliptic_hecke_jets_stage2_data.m').read_text();ss=lit(b,'SourceRecords');tt=lit(b,'SourceHeckeRecords');F=GF(2**18,'z');R=PolynomialRing(F,'X');X=R.gen();FB=Field(bits([int(c)for c in F.modulus()]));assert FB.d==18
code=lambda a:bits([int(c)for c in a.polynomial()])
def units(E,P):
 a1,a2,a3,a4,a6=E.a_invariants();x,y=P[0],P[1];aa3=a1*x+a3;aa2=x+a2;aa4=x*x+a4+a1*y;aaa1=a1;w=aa4/aa3;aaa2=aa2+w*w+a1*w;b=aaa2**3/aa3**2;c=(aaa1*aaa2+aa3)/aa3;r=b/c;s=c*c/(b+c);v=r*s+1;return((s+r)/v,v/(s*s+s+r+1))
results=[]
for sid,(s,trs)in enumerate(zip(ss,tt),1):
 v=decode_record(s);roots=(R(v['field_poly'])).roots(multiplicities=False);assert len(roots)==9;r=roots[0]
 def embed(a):return sum(F((a>>j)&1)*r**j for j in range(9))
 E=EllipticCurve(F,list(map(embed,v['ainv'])));P=E(list(map(embed,v['point'])));assert P.order()==85;pol=X**4+X**3+E.a6();ks=pol.roots(multiplicities=False);assert len(ks)==4;actual=collections.Counter()
 for x in ks:
  ys=(X**2+x*X-(x**3+E.a2()*x*x+E.a4()*x+E.a6())).roots(multiplicities=False);assert len(ys)==2;Q=E([x,ys[0]]);assert Q.order()==3;phi=E.isogeny([Q]);assert phi.degree()==3;ET=phi.codomain();PT=phi(P);assert PT.order()==85;xy=units(ET,PT);d=9 if all(FB.frob(code(a),9)==code(a)for a in xy)else 18;k=FB.key([1]+[code(a)for a in xy],d);actual[k]+=1
 expected=collections.Counter()
 for e,d,m,rec in trs:
  w=decode_record(rec);expected[w['F'].key([1]+w['units'],d)]+=e
 assert actual==expected,(sid,actual,expected);results.append({'source':sid,'kernel_polynomials':4,'geometric_targets':4,'closed_target_degrees':[9,9,18],'level_order':85});print('NONIC85_SAGE_VELU_TRANSPORT_SOURCE_PASS',sid,flush=True)
out={'status':'PASS','sources':results,'all_twenty_level_transports_independently_rebuilt':True}
p=H/'hecke_transport_independent_verified.json'
if '--readonly'in sys.argv:assert json.loads(p.read_text())==out
else:p.write_text(json.dumps(out,indent=2)+'\n')
print('NONIC85_SAGE_VELU_TRANSPORT_PASS elapsed',time.monotonic()-t)
