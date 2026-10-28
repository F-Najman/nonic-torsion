#!/usr/bin/env python3
"""Independent finite-field and integer-lattice replay; no Magma is invoked."""
from pathlib import Path
from collections import Counter
import json,math,sys,time
from sage.all import GF,PolynomialRing,EllipticCurve,matrix,ZZ,identity_matrix
P=Path(__file__).resolve().parent.parent/'data'/'C62_X131'

def closed_counts(pc,q,n=9):
 a=list(reversed(pc));assert a[0]==1
 power=[0];closed=[]
 for k in range(1,n+1):
  power.append(-k*a[k]-sum(a[j]*power[k-j] for j in range(1,k)))
  v=q**k+1-power[k]-sum(d*closed[d-1] for d in range(1,k) if k%d==0)
  assert v%k==0;closed.append(v//k)
 return closed

def main():
 started=time.monotonic();W=json.loads((P/'presentation/presentation.json').read_text());D=json.loads((P/'characters62.json').read_text());C=json.loads((P/'classification62.json').read_text());G=json.loads((P/'geometry62.json').read_text())
 assert W['native_proven'] and W['checkpoint_complete'] and W['factor_basis_bound']==9
 assert (W['N'],W['p'],W['genus'],W['field_degree'])==(31,3,26,13)
 assert D['complete'] and C['complete'] and G['complete']
 m=3481;fb=W['factor_basis'];assert len(fb)==m and [r['index'] for r in fb]==list(range(1,m+1))
 deg=[r['degree'] for r in fb];expected=[15,0,0,15,91,120,330,930,1980]
 assert [deg.count(i) for i in range(1,10)]==expected==closed_counts(W['frobenius_coefficients'],3)
 assert W['frobenius_coefficients']==G['N31_frobenius_coefficients'] and sum(W['frobenius_coefficients'])==W['known_Pic_order']==134097116106800
 assert W['candidate_invariants']==[2,2,10,3352427902670,0] and math.prod(W['candidate_invariants'][:-1])==W['known_Pic_order']
 for co,g in ((W['frobenius_coefficients'],26),(G['N62_frobenius_coefficients'],91)):
  assert len(co)==2*g+1 and all(co[i]==3**(g-i)*co[2*g-i] for i in range(g+1))
 source=closed_counts(G['N62_frobenius_coefficients'],3);source[0]-=30
 target=expected[:];target[0]-=15
 assert source==G['N62_noncuspidal_counts'] and source[:4]==[0]*4
 assert target==G['N31_noncuspidal_counts']==C['noncuspidal_counts']==[0,0,0,15,91,120,330,930,1980]
 base=D['base_FB_index']-1;cusps=D['cusp_FB_indices'];assert len(cusps)==len(set(cusps))==15 and sorted(cusps)==[i+1 for i,d in enumerate(deg) if d==1] and base+1 in cusps
 Bs=[D['matrix16'],D['matrix31'],D['matrix61']];mods=[16,31,61];full=[16,31,61];horders=[4,1,1]
 assert all(len(row)==m and row[base]==0 for B in Bs for row in B)
 for rel in W['relations']:
  assert len(rel)==m;sparse=[(i,int(x)) for i,x in enumerate(rel) if x]
  assert sum(deg[i]*x for i,x in sparse)==0
  for B,mod in zip(Bs,mods):assert all(sum(x*b[i] for i,x in sparse)%mod==0 for b in B)
 assert len(W['relations'])==4391
 def lattice(B,mod,indices=None):
  A=matrix(ZZ,B);V=A if indices is None else A.matrix_from_columns(indices)
  return V.augment(mod*identity_matrix(ZZ,A.nrows())).transpose().hermite_form(include_zero_rows=False).transpose()
 tests=[]
 for B,mod,fo,ho in zip(Bs,mods,full,horders):
  r=len(B);I=lattice(B,mod);H=lattice(B,mod,[i-1 for i in cusps])
  assert I.nrows()==I.ncols()==r and mod**r//abs(I.det())==fo
  assert H.nrows()==H.ncols()==r and mod**r//abs(H.det())==ho
  A=mod*H.inverse();assert all(x in ZZ for x in A.list());tests.append(matrix(ZZ,A))
 assert math.prod(full)//math.prod(horders)==7564 and W['known_Pic_order']//7564==10*1772833370
 print('C62_RELATIONS_AND_PRIMARY_IMAGES_PASS',len(W['relations']),flush=True)
 k=GF(3);R=PolynomialRing(k,'z');fields={};common={d:GF(3**d,name='c') for d in range(3,10)};embeddings={};points={};keys=set();count=Counter();lifts=Counter();eligible=Counter();rootcounts=Counter()
 lines=(P/'finite_points62.jsonseq').read_text().splitlines();header=json.loads(lines[0]);finish=json.loads(lines[-1]);assert header['FB_size']==m and finish['type']=='complete'
 def evalpoly(coeff,a):
  v=a.parent()(0)
  for c in reversed(coeff):v=v*a+int(c)
  return v
 for line in lines[1:-1]:
  a=json.loads(line);assert a['type']=='point';i=a['FB_index'];d=a['degree'];assert i not in points and deg[i-1]==d and d>=3
  cp=tuple(a['field_polynomial'])
  if cp not in fields:
   f=R(list(cp));assert f.is_irreducible() and f.degree()==d
   fields[cp]=GF(3**d,name='v',modulus=f)
   roots=f.change_ring(common[d]).roots();assert len(roots)==d
   embeddings[cp]=fields[cp].hom([roots[0][0]],common[d])
  L=fields[cp];v=L.gen();b=evalpoly(a['b'],v);c=evalpoly(a['c'],v)
  assert len(a['b'])==len(a['c'])==d and b and c and b!=c
  E=EllipticCurve(L,[1-c,-b,-b,0,0]);Q=E(0,0);assert not E.is_singular() and 31*Q==E(0) and Q!=E(0)
  S=PolynomialRing(L,'Z');z=S.gen();a1,a2,a3,a4,a6=E.ainvs();two=4*z**3+(a1*a1+4*a2)*z*z+2*(a1*a3+2*a4)*z+a3*a3+4*a6
  assert two.gcd(two.derivative()).degree()==0
  nr=len(two.roots());assert nr==a['rational_two_points'] and nr in (0,1,3)
  if d in (3,4):assert nr==0 and two.is_irreducible()
  # Invert the Tate normal form coordinates, then check the MD affine equation.
  r=b/c;s=c*c/(b-c);T=r*s-2*r+1;x=(s-r)/T;y=T/(s*s-s-r+1)
  rr=(x*x*y-x*y+y-1)/(x*x*y-x);ss=(x*y-y+1)/(x*y)
  assert (rr*ss*(rr-1),ss*(rr-1))==(b,c)
  assert sum(int(co)*x**e[0]*y**e[1] for e,co in zip(G['N31_affine_exponents'],G['N31_affine_coefficients']))==0
  t=y;u=x*(2*t**4+t**7);cache={}
  def poly(c):
   key=tuple(c)
   if key not in cache:cache[key]=evalpoly(c,t)
   return cache[key]
  def rat(c):
   den=poly(c[1]);assert den
   return poly(c[0])/den
  def elt(c):
   w=L(0)
   for cr in reversed(c):w=w*u+rat(cr)
   return w
  assert sum(rat(co)*u**j for j,co in enumerate(W['field_polynomial']))==0
  assert fb[i-1]['finite'] and len(fb[i-1]['ideal_basis'])==13
  assert all(elt(row)==0 for row in fb[i-1]['ideal_basis'])
  # Compare in one common field per degree, allowing different primitive polynomials.
  orb=[];xx=embeddings[cp](x);yy=embeddings[cp](y);xcommon,ycommon=xx,yy
  for j in range(d):
   orb.append((tuple(int(e) for e in xx.polynomial().list()),tuple(int(e) for e in yy.polynomial().list())))
   xx=xx**3;yy=yy**3
  assert (xx,yy)==(xcommon,ycommon) and len(set(orb))==d
  key=(d,min(orb));assert key not in keys;keys.add(key)
  count[d]+=1;lifts[d]+=nr;rootcounts[d,nr]+=1
  mult=3 if d==3 else (1 if d>=5 and nr else 0)
  if mult:eligible[d]+=1
  points[i]=(d,mult);a['computed_mult']=mult
  if len(points)%250==0:print('C62_POINT_REPLAY',len(points),'seconds',round(time.monotonic()-started,2),flush=True)
 assert [count[i] for i in range(1,10)]==target and set(points)|set(cusps)==set(range(1,m+1))
 covered=[lifts[i] for i in range(1,10)];covered[8]+=count[3]
 assert covered==source
 expectedtests={i:(d,mult) for i,(d,mult) in points.items() if mult}
 got={};survivors=[]
 for T in C['tests']:
  i=T['FB_index'];assert i not in got;got[i]=True;d,mult=expectedtests[i]
  assert (T['degree'],T['multiplicity'])==(d,mult)
  images=[[mult*row[i-1]%mod for row in B] for B,mod in zip(Bs,mods)]
  assert images==[T['image16'],T['image31'],T['image61']]
  member=all(all(sum(A[j,h]*v[h] for h in range(A.ncols()))%mod==0 for j in range(A.nrows())) for A,v,mod in zip(tests,images,mods))
  assert member==T['in_rational_image']
  if member:survivors.append(i)
 assert set(got)==set(expectedtests) and sorted(survivors)==sorted(C['survivor_indices'])
 assert [eligible[i] for i in range(1,10)]==C['eligible_counts']==finish['eligible_counts']
 assert finish['tests']==len(got) and finish['survivors']==len(survivors) and finish['noncuspidal_counts']==target
 report={'level':31,'candidate_cyclic_order':62,'degree':9,'relations':4391,'FB_places':3481,'finite_points_verified':3466,'field_and_all_ideal_basis_bindings':True,'distinct_Frobenius_orbits':3466,'complete_target_and_source_counts':True,'primary_image_orders':full,'rational_primary_image_orders':horders,'rational_image_index':7564,'eligible_counts':C['eligible_counts'],'tests':len(got),'survivor_FB_indices':sorted(survivors),'class_sieve_verified':True,'complete':True,'elapsed_seconds':round(time.monotonic()-started,3)}
 print('NONIC62_INDEPENDENT_REPLAY '+json.dumps(report),flush=True)
 if '--write-report' in sys.argv:(P/'independent_replay62.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
