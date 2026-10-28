"""Independent finite-field point counts and classical Phi_5 replay."""
from sage.all import *
from sage.schemes.elliptic_curves.mod_poly import classical_modular_polynomial
from pathlib import Path
import json,re,ast,sys
p=Path(__file__).resolve().parents[1]/'data'/'direct5'
s=(p/'isogeny5_target_data.m').read_text();data={}
for name,value in re.findall(r'Append\(~(\w+),(.*?)\);',s,re.S):data.setdefault(name,[]).append(ast.literal_eval(value))
assert len(data['SourceDegrees'])==4 and data['SourceDegrees']==[9]*4
assert data['TargetDegrees']==[54]*4
assert data['SourceTargetIndices']==[[1],[2],[3],[4]]
# Direct ordinary point count: for y^2+xy=x^3+a2*x^2+1/j,
# each x!=0 has 1+(-1)^Tr(x+a2+1/(j*x^2)) points.
census=[]
for d in [8,9]:
 F=GF(2**d,'v');R=PolynomialRing(F,'T');orbits=set();qual=0
 for j in F:
  if not j:continue
  sm=sum(1 if (xx+1/(j*xx*xx)).trace()==0 else -1 for xx in F if xx)
  counts=[int(2**d+1+sm),int(2**d+1-sm)]
  for n in counts:
   if n%87==0:
    qual+=1;mp=j.minpoly();assert mp.degree()==d
    orbits.add(tuple(int(v) for v in mp.list()))
 assert len(orbits)==(0 if d==8 else 4)
 assert qual==d*len(orbits)
 if d==9:assert orbits=={tuple(r) for r in data['SourceJPolynomials']}
 # Every j=0 curve is isomorphic to y^2+A*y=x^3+B*x+C:
 # translate x to remove x^2, scale A through cube classes, and
 # translate y to reduce C/A^2 to its two Artin--Schreier classes.
 supersingular_orders=set();g=F.multiplicative_generator()
 for k in range(gcd(3,2**d-1)):
  A=g**k
  for B in F:
   sm0=sum(1 if ((xx**3+B*xx)/A**2).trace()==0 else -1 for xx in F)
   supersingular_orders.update([int(2**d+1+sm0),int(2**d+1-sm0)])
 assert all(n%87 for n in supersingular_orders)
 census.append({'degree':d,'ordinary_j_checked':2**d-1,'qualifying_j_twists':qual,'qualifying_j_orbits':len(orbits),'jzero_orders':sorted(supersingular_orders),'jzero_all_normal_forms_checked':True})
 print('DIRECT_ORDINARY_CENSUS',census[-1],flush=True)
F2=GF(2);R2=PolynomialRing(F2,'T')
M=GF(2**54,'w');RM=PolynomialRing(M,'Z');Z=RM.gen()
phi=classical_modular_polynomial(5)
odd=sorted((tuple(map(int,e)),int(v)%2) for e,v in phi.dict().items() if int(v)%2)
assert [e for e,v in odd]==[(0,6),(2,4),(4,2),(5,5),(6,0)]
records=[]
for i in range(4):
 L=GF(2**9,'a',modulus=R2(data['SourceFieldPolynomials'][i]));a=L.gen()
 ev=lambda co:sum(L(v)*a**k for k,v in enumerate(co))
 E=EllipticCurve(L,[ev(co) for co in data['SourceECoefficients'][i]])
 P=E([ev(co) for co in data['SourcePointCoefficients'][i]])
 assert E.discriminant()!=0 and 87*P==E(0) and 29*P!=E(0) and 3*P!=E(0)
 assert [int(v) for v in E.j_invariant().minpoly().list()]==data['SourceJPolynomials'][i]
 # Independent literal point count on the saved curve.
 A1,A2,A3,A4,A6=E.a_invariants();count=1
 for x in L:
  B=A1*x+A3;C=x**3+A2*x*x+A4*x+A6
  count+=1 if B==0 else (2 if (C/B**2).trace()==0 else 0)
 assert count==522
 sr=RM(data['SourceJPolynomials'][i]).roots();assert len(sr)==9
 j=sr[0][0];F=Z**6+j**5*Z**5+j**2*Z**4+j**4*Z**2+j**6
 rr=F.roots();assert len(rr)==6 and all(mult==1 for r,mult in rr)
 target=R2(data['TargetJPolynomials'][data['SourceTargetIndices'][i][0]-1]);assert target.degree()==54 and target.is_irreducible()
 assert all(R2(r.minpoly())==target for r,mult in rr)
 # The six roots are precisely one F_(2^9)-Frobenius orbit.
 z=rr[0][0];assert len({z**(2**(9*k)) for k in range(6)})==6
 assert {z**(2**(9*k)) for k in range(6)}=={r for r,mult in rr}
 records.append({'source_index':i+1,'source_degree':9,'target_degree':54,'six_distinct_isogeny_j_values':True,'Frob9_cycle_length':6,'norm_multiplicity':1,'point_order':87,'curve_order':count})
 print('PHI5_AND_SOURCE_POINT',records[-1],flush=True)
out={'complete':True,'ordinary_census':census,'phi5_mod2_monomials':[[int(a),int(b)] for (a,b),v in odd],'records':records}
if '--check' in sys.argv:assert json.loads((p/'independent_replay5.json').read_text())==out
else:(p/'independent_replay5.json').write_text(json.dumps(out,indent=2)+'\n')
print('NONIC87_INDEPENDENT_PHI5_REPLAY_COMPLETE',flush=True)
