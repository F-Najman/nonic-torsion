from pathlib import Path
import ast,re,json,sys,time
from sage.all import GF,PolynomialRing,EllipticCurve,gcd
P=Path(__file__).resolve().parents[1]/'data'/'CM6';R=P.parents[2];t0=time.monotonic()
def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']'))
def named(s,n):return lit(re.search(r'\b'+n+r'\s*:=\s*(\[.*?\]);',s,re.S)[1])
def rd(f):return lit((P/f).read_text())
pr=PolynomialRing(GF(2),'a');out=[];native=rd('good6_jets_data.txt')[2][0][0]
for d in range(6,10):
 s=(R/'C65'/'data'/'census'/f'N65_f{d}_records.m').read_text();pol=named(s,'FieldPolynomial');K=GF(2**d,'a',modulus=pr(pol),impl='pari_ffelt');a=K.gen();q=K.cardinality();values=list(K);codes={v:int(sum(int(c)<<j for j,c in enumerate(v._vector_())))for v in values};values.sort(key=codes.get);gg=K.multiplicative_generator();cubeclasses=gcd(3,q-1);cases=[];eligible=[]
 for u in range(cubeclasses):
  a3=gg**u;image={v**4+a3*v for v in values};todo=set(values);ar=[]
  while todo:
   v=min(todo,key=codes.get);ar.append(v);todo-={v+w for w in image}
  assert len(ar)*len(image)==q
  im6={v*v+a3*v for v in values};nonzero=next(v for v in values if v not in im6);assert len(im6)==q//2
  for a4 in ar:
   for a6 in[K(0),nonzero]:
    E=EllipticCurve(K,[0,0,a3,a4,a6]);order=1+2*sum(((x**3+a4*x+a6)/a3**2).trace()==0 for x in values);assert order==E.cardinality();cases.append(order)
    if order%65==0:eligible.append(E)
 if d==6:
  assert eligible
  def dec(c):return sum(K(v)*a**j for j,v in enumerate(c))
  assert native[0]==pol;E0=EllipticCurve(K,[dec(c)for c in native[1]]);PP=E0(*[dec(c)for c in native[2]]);assert E0.cardinality()==65 and 65*PP==E0(0)and 13*PP!=E0(0)and 5*PP!=E0(0);assert all(E.is_isomorphic(E0)for E in eligible)
  assert E0.abelian_group().invariants()==(65,)
  eligible_orders=sorted(set(int(E.cardinality())for E in eligible))
 else:assert not eligible;eligible_orders=[]
 out.append(dict(degree=d,supersingular_containing_models=len(cases),qualified_models=len(eligible),qualified_orders=eligible_orders,all_j0_coefficients_covered_by_cube_additive_cosets=True))
 print('C65_SUPERSINGULAR_CENSUS_DEGREE_PASS',d,len(cases),len(eligible),flush=True)
report=dict(complete=True,level=65,prime=2,all_supersingular_degrees_6_through_9=out,unique_eligible_F64_isomorphism_class_matches_native=True,eligible_F64_group_cyclic65=True)
if '--readonly'in sys.argv:assert json.loads((P/'supersingular_census_sage_certificate.json').read_text())==report
else:(P/'supersingular_census_sage_certificate.json').write_text(json.dumps(report,indent=2)+'\n')
print('C65_SUPERSINGULAR_COMPLETE_CENSUS_INDEPENDENT_PASS',json.dumps(report),'seconds',time.monotonic()-t0,flush=True)
