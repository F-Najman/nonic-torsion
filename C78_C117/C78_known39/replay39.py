from sage.all import *
from pathlib import Path
import json,collections,sys,hashlib,re
p=Path(__file__).resolve().parents[1]/'data/C78_known39';J=lambda f:json.loads((p/f).read_text());z=J('inputs/known39.json')['example'];g=J('global39.json');raw=(p/'diamond39.jsonseq').read_text();dec=json.JSONDecoder();d=[]
while raw.strip():raw=raw.lstrip();v,i=dec.raw_decode(raw);d.append(v);raw=raw[i:]
R=PolynomialRing(QQ,'T');T=R.gen();f=R([QQ(a)for a in z['f']]);xp=R([QQ(a)for a in z['x']]);yp=R([QQ(a)for a in z['y']]);assert f.degree()==9 and f.is_irreducible()and f.discriminant()==-3**14*31**3
nativeinput=(p/'inputs/known39.m').read_text()
for field,key in [('FieldCoefficients','f'),('XCoefficients','x'),('YCoefficients','y')]:
 terms=re.search(field+r':=\[(.*?)\];',nativeinput,re.S)[1].split(',');assert [QQ(c)for c in terms]==[QQ(c)for c in z[key]]
K=NumberField(f,'a');a=K.gen();x=xp(a);y=yp(a)
def fromxy(x,y):
 r=(x*x*y-x*y+y-1)/(x*x*y-x);s=(x*y-y+1)/(x*y);b=r*s*(r-1);c=s*(r-1);return EllipticCurve([1-c,-b,-b,0,0])
def tatebc(E,P):
 r,t=P[0],P[1];a1,a2,a3,a4,a6=E.a_invariants();a3p=a3+2*t+a1*r;s=(3*r*r+2*a2*r+a4-a1*t)/a3p
 A1=a1+2*s;A2=a2+3*r-s*a1-s*s;A3=a3p;u=A3/A2;b=-A2/u**2;c=1-A1/u
 A4=a4+2*r*a2+3*r*r-s*a3-(t+r*s)*a1-2*s*t;A6=a6+r*a4+r*r*a2+r**3-t*a3-r*t*a1-t*t
 assert A1/u==1-c and A2/u**2==-b and A3/u**3==-b and A4==0 and A6==0
 return b,c
def qqpair(v):return QQ(v[0])/QQ(v[1])
def nf(v):return sum(qqpair(c)*a**i for i,c in enumerate(v))
E=fromxy(x,y);P=E(0,0);assert 39*P==E(0)and 13*P!=E(0)and 3*P!=E(0)and E.j_invariant().minpoly().degree()==3;assert list(E.a_invariants())==[nf(v)for v in g['global_A_invariants']]
assert [qqpair(v)for v in g['field_polynomial']]==f.list()and g['polynomial_discriminant']==f.discriminant()
assert g['point_order']==39 and g['field_degree']==9 and g['global_trace_groups']and len(d)==12
units=[i for i in range(1,20)if gcd(i,39)==1];assert units==g['physical_multipliers']==[v['multiplier']for v in d]
k=GF(2);R2=PolynomialRing(k,'t');f2=R2(f);assert f2==R2([1]+[0]*7+[1,1])and f2.is_irreducible();L=k.extension(f2,'z');b0=L.gen();vec=lambda v:[int(L(v).polynomial()[i])for i in range(9)]
def reduce_nf(v,field,root):return sum((field(c.numerator())/field(c.denominator()))*root**i for i,c in enumerate(K(v).list()))
def key(v):return min(tuple(vec(v[0]**(2**i))+vec(v[1]**(2**i)))for i in range(9))
x2=R2(xp)(b0);y2=R2(yp)(b0);E2=fromxy(x2,y2);P2=E2(0,0);assert E2.discriminant()!=0 and 39*P2==E2(0)and 13*P2!=E2(0)and 3*P2!=E2(0)
count2=1;aa=E2.a_invariants()
for xx in L:
 A=aa[0]*xx+aa[2];B=xx**3+aa[1]*xx*xx+aa[3]*xx+aa[4]
 count2+=(1 if A==0 else 2 if (B/(A*A)).trace()==0 else 0)
assert count2==g['finite2_order']==468
roots=[nf(v)for v in g['global_automorphism_generators']];assert len(set(roots))==3 and all(f(v)==0 for v in roots)
# Every automorphism group has order dividing9. Negative discriminant rulesout a Galois degree9 field, so these three automorphisms are all of them.
assert f.discriminant()<0
apply=lambda r,v:sum(c*r**i for i,c in enumerate(K(v).list()))
FKeys=[];GBC=[]
for v in d:
 m=v['multiplier'];gb=tatebc(E,m*P);assert list(gb)==[nf(c)for c in v['global_Tate_bc']];assert gb==tatebc(E,(39-m)*P)
 lam=v['primitive_lambda'];h=R([qqpair(q)for q in v['primitive_polynomial']]);assert h.degree()==9 and h.is_irreducible()and h(gb[0]+lam*gb[1])==0
 fb=tatebc(E2,m*P2);assert len({tuple(q**(2**i)for q in fb)for i in range(9)})==9;assert [vec(q)for q in fb]==v['finite2_Tate_bc']and [reduce_nf(q,L,b0)for q in gb]==list(fb)
 kk=key(fb);assert list(kk)==v['finite2_frobenius_key'];FKeys.append(kk);GBC.append(gb)
for i,v in enumerate(d):
 orb=[units[j]for j,bb in enumerate(GBC)if any(tuple(apply(r,q)for q in GBC[i])==bb for r in roots)];finite=[units[j]for j,kk in enumerate(FKeys)if kk==FKeys[i]];assert orb==finite==v['trace_group']and len(orb)==3
assert len(set(FKeys))==4 and {tuple(v)for v in g['global_trace_groups']}=={tuple(v['trace_group'])for v in d}
exponents=[]
for rt in roots:
 red=reduce_nf(rt,L,b0);ee=[i for i in range(9)if red==b0**(2**i)];assert len(ee)==1;exponents.append(ee[0])
assert exponents==g['automorphism_Frobenius_exponents']and sorted(exponents)==[0,3,6]
p29=GF(29);x29=p29(10);R29=PolynomialRing(p29,'U');assert R29(f)(x29)==0
def eval29(poly):return sum(p29(c.numerator())/p29(c.denominator())*x29**i for i,c in enumerate(poly.list()))
x29=eval29(xp);y29=eval29(yp);assert(x29,y29)==(10,20);E29=fromxy(x29,y29);assert list(map(int,E29.a_invariants()))==[3,12,12,0,0]and E29.discriminant()==1
A1,A2,A3,A4,A6=E29.a_invariants();count29=1+sum(1 for xx in p29 for yy in p29 if yy*yy+A1*xx*yy+A3*yy==xx**3+A2*xx*xx+A4*xx+A6);assert count29==39 and f.discriminant()%29!=0
Q=E29(0,0);assert 39*Q==E29(0)and 13*Q!=E29(0)and 3*Q!=E29(0)
known=[]
for kk in sorted(set(FKeys)):
 members=[v for v in d if tuple(v['finite2_frobenius_key'])==kk];v=min(members,key=lambda v:v['multiplier'])
 known.append({'representative_multiplier':v['multiplier'],'physical_multipliers':v['trace_group'],'global_trace_degree':9,'reduction_code':[9,list(map(int,f2.list())),v['finite2_Tate_bc']],'frobenius_canonical_Tate_bc_key':list(kk)})
binding={'complete':True,'source_example_index':16,'prime':2,'original_F512_polynomial':list(map(int,f2.list())),'global_closed_nonic_divisors':4,'physical_multipliers':12,'each_finite_orbit_multiplicity':3,'known_positive_divisors':known}
out={'complete':True,'global_exact_point_order':39,'full_global_torsion_order':39,'good_unramified_prime':29,'rational_residue_root':10,'F29_direct_pairs_enumerated':841,'F29_order':39,'F512_trace_count':468,'global_trace_orbits':4,'physical_multiplier_bindings':12,'all_global_automorphisms_checked':3,'all_Frobenius_keys_checked':True}
if'--check'in sys.argv:
 assert J('positive_divisors39.json')==binding and J('independent39.json')==out
else:
 (p/'positive_divisors39.json').write_text(json.dumps(binding,indent=2)+'\n');(p/'independent39.json').write_text(json.dumps(out,indent=2)+'\n')
print('NONIC78_KNOWN39_INDEPENDENT_REPLAY_PASS exact39_full39_F29_direct=true global4_finite4_bindings=true',flush=True)
