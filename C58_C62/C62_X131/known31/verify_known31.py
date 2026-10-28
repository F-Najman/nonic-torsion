"""Exact positive degree-nine binding.  Run with sage -python.
This script does not assume the finite class-filter conclusion; it binds its
five stated survivors to five reductions of actual global diamond points.
"""
from sage.all import *
from pathlib import Path
import hashlib,json,sys,time
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
TOP=ROOT/'C58_C62'/'data'/'C62_X131'
INP=ROOT/'common'/'data'/'inputs'/'vanhoeij_nonic_input.json'
assert INP.is_file(),str(INP)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(INP)=='fed7f3544fbea4b88cebdfb9996bec5549de4a499dfefc81960378254a1f8ce5'
ex=next(e for e in json.loads(INP.read_text())['examples'] if e['index']==2)
assert ex['N']==31
RQ=PolynomialRing(QQ,'X'); X=RQ.gen()
f=RQ([QQ(s) for s in ex['f']]);assert f.is_irreducible()
K=NumberField(f,'a');a=K.gen()
x=sum(QQ(s)*a**i for i,s in enumerate(ex['x']))
y=sum(QQ(s)*a**i for i,s in enumerate(ex['y']))
def bc_xy(x,y):
 r=(x*x*y-x*y+y-1)/(x*x*y-x);s=(x*y-y+1)/(x*y)
 return r*s*(r-1),s*(r-1)
def xy_bc(b,c):
 r=b/c;s=c*c/(b-c);T=r*s-2*r+1
 return (s-r)/T,T/(s*s-s-r+1)
def tate(E,P):
 r,t=P[0],P[1];a1,a2,a3,a4,a6=E.a_invariants()
 ap=a3+2*t+a1*r
 m=(3*r*r+2*a2*r+a4-a1*t)/ap
 A1=a1+2*m;A2=a2+3*r-m*a1-m*m;A3=ap;u=A3/A2
 return -A2/u**2,1-A1/u
b,c=bc_xy(x,y);assert xy_bc(b,c)==(x,y)
E=EllipticCurve(K,[1-c,-b,-b,0,0]);P=E(0,0)
assert E.discriminant()!=0 and P!=E(0) and 31*P==E(0)
assert tate(E,P)==(b,c)
# Good degree-one reductions at distinct odd primes force full torsion 31.
reductions=[]
for p,target in [(53,62),(83,93)]:
 F=GF(p);pol=PolynomialRing(F,'z')([F(v) for v in f.list()]);found=None
 for ar,mul in pol.roots():
  if mul!=1:continue
  try:
   def ev(coeff):
    return sum(F(QQ(s).numerator())/F(QQ(s).denominator())*ar**i for i,s in enumerate(coeff))
   xr,yr=ev(ex['x']),ev(ex['y']);br,cr=bc_xy(xr,yr)
   Er=EllipticCurve(F,[1-cr,-br,-br,0,0])
   if Er.discriminant()==0:continue
   Pr=Er(0,0);assert Pr!=Er(0) and 31*Pr==Er(0)
   aa=Er.a_invariants(); n=1
   for xx in F:
    disc=(aa[0]*xx+aa[2])**2+4*(xx**3+aa[1]*xx**2+aa[3]*xx+aa[4])
    n+=1 if disc==0 else (2 if disc.is_square() else 0)
   if n==target:found={'p':p,'root':int(ar),'root_simple':True,'good':True,'order':int(n)};break
  except ZeroDivisionError:continue
 assert found is not None,(p,target)
 reductions.append(found)
assert gcd([r['order'] for r in reductions])==31
print('GLOBAL31_EXACT_ORDER_AND_FULL_TORSION',reductions,flush=True)
F3=GF(3);R3=PolynomialRing(F3,'z');fm=R3([F3(v) for v in f.list()]);assert fm.is_irreducible()
L=GF(3**9,'z',modulus=fm);z=L.gen()
def enc(v):return tuple(int(q) for q in v.polynomial().list())+(0,)*(9-len(v.polynomial().list()))
def evQ(coeff):return sum(L(QQ(s).numerator())/L(QQ(s).denominator())*z**i for i,s in enumerate(coeff))
x3,y3=evQ(ex['x']),evQ(ex['y']);b3,c3=bc_xy(x3,y3)
E3=EllipticCurve(L,[1-c3,-b3,-b3,0,0]);P3=E3(0,0)
assert E3.discriminant()!=0 and 31*P3==E3(0) and P3!=E3(0)
RR=PolynomialRing(L,'w');w=RR.gen();aa=E3.a_invariants()
twopoly=4*(w**3+aa[1]*w**2+aa[3]*w+aa[4])+(aa[0]*w+aa[2])**2
assert len(twopoly.roots())==3
geom=json.loads((TOP/'geometry62.json').read_text())
def curve31(xx,yy):
 return sum(L(v)*xx**e[0]*yy**e[1] for e,v in zip(geom['N31_affine_exponents'],geom['N31_affine_coefficients']))
assert curve31(x3,y3)==0
# Frobenius orbit of ordered Tate coefficients determines a closed fine-moduli point.
def key(bb,cc):
 return min((enc(bb**(3**i)),enc(cc**(3**i))) for i in range(9))
def flat(k):return list(k[0]+k[1])
known={};multipliers=[]
for k in range(1,16):
 bb,cc=tate(E3,k*P3);xx,yy=xy_bc(bb,cc)
 assert curve31(xx,yy)==0
 assert all((bb**(3**j),cc**(3**j))!=(bb,cc) for j in [1,3])
 EE=EllipticCurve(L,[1-cc,-bb,-bb,0,0]);PP=EE(0,0)
 assert 31*PP==EE(0) and PP!=EE(0) and EE.j_invariant()==E3.j_invariant()
 kk=key(bb,cc);known.setdefault(kk,[]).append(k)
 multipliers.append({'multiplier':k,'b':list(enc(bb)),'c':list(enc(cc)),'orbit_key':flat(kk)})
assert len(known)==5 and sorted(map(len,known.values()))==[3]*5
cl=json.loads((TOP/'classification62.json').read_text());assert cl['complete']
surv=cl['survivor_indices'];assert surv==[1878,2248,3105,3225,3323]
points=[json.loads(l) for l in (TOP/'finite_points62.jsonseq').read_text().splitlines()]
points={v['FB_index']:v for v in points if v.get('type')=='point' and v.get('FB_index') in surv}
bindings=[];seen=set()
for idx in surv:
 rec=points[idx];pol=RR([L(v) for v in rec['field_polynomial']]);roots=sorted([r for r,m in pol.roots()],key=enc)
 assert len(roots)==9
 root=roots[0];bb=sum(L(v)*root**i for i,v in enumerate(rec['b']));cc=sum(L(v)*root**i for i,v in enumerate(rec['c']))
 xx,yy=xy_bc(bb,cc);assert curve31(xx,yy)==0
 kk=key(bb,cc);assert kk in known and kk not in seen;seen.add(kk)
 EE=EllipticCurve(L,[1-cc,-bb,-bb,0,0]);assert EE.j_invariant() in [E3.j_invariant()**(3**i) for i in range(9)]
 assert 31*EE(0,0)==EE(0)
 bindings.append({'FB_index':idx,'degree':9,'embedding_root':list(enc(root)),'orbit_key':flat(kk),'known_global_diamond_multipliers':known[kk]})
assert seen==set(known)
out={'complete':True,'known_example_index':2,'level':31,'field_degree':9,'field_polynomial':ex['f'],'input_sha256':sha(INP),'classification_sha256':sha(TOP/'classification62.json'),'finite_points_sha256':sha(TOP/'finite_points62.jsonseq'),'geometry_sha256':sha(TOP/'geometry62.json'),'good_reductions':reductions,'full_torsion':'C31','reduction_mod3_irreducible':True,'residue_two_roots':3,'common_finite_field_polynomial':[int(v) for v in fm.list()],'multipliers':multipliers,'survivor_bindings':bindings}
if '--check' in sys.argv:
 assert json.loads((TOP/'known31'/'binding31.json').read_text())==out
else:(TOP/'known31'/'binding31.json').write_text(json.dumps(out,indent=2)+'\n')
print('KNOWN31_FIVE_SURVIVORS_BOUND',[(r['FB_index'],r['known_global_diamond_multipliers']) for r in bindings],flush=True)
print('KNOWN31_BINDING_COMPLETE full_torsion=31 survivor_count=5',flush=True)
