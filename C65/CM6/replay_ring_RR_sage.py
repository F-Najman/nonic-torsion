from pathlib import Path
import ast,re,json,time,sys,hashlib
from sage.all import GF,Matrix,vector,PolynomialRing
P=Path(__file__).resolve().parents[1]/'data'/'CM6';t0=time.monotonic();F=GF(2);g=121;h=360;prec=1100

def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']'))
def rd(f):return lit((P/f).read_text())
def mat(s,name):
 m=re.search(name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,name
 return Matrix(F,int(m[1]),int(m[2]),lit(m[3]))
def unhex(s,n):
 k=int(s,16);assert k>>n==0;return[(k>>i)&1 for i in range(n)]
source=(P.parent/'source_units'/'canonical_unit_data.m').read_bytes();sha=hashlib.sha256(source).hexdigest();Q=mat(source.decode(),'qcoeffs');assert Q.dimensions()==(g,prec)and Q[:,:337].rank()==g
T=json.loads((P.parent/'composite65'/'canonical_ring_R1_tensor.json').read_text());assert T['native_data_sha256']==sha and T['dimensions']==[g,h]and T['q_precision']==prec and T['Sturm_weight4']==672
Q2=Matrix(F,[unhex(x,prec)for x in T['basis_q']]);assert Q2.dimensions()==(h,prec)and Q2[:,:673].rank()==h
pairs=[(i,j)for i in range(g)for j in range(i,g)];X=T['product_1_1_into_2'];assert len(X)==g and all(len(row)==g for row in X)and all(X[i][j]==X[j][i]for i in range(g)for j in range(g))
Tensor=Matrix(F,[unhex(X[i][j],h)for i,j in pairs]);assert Tensor.dimensions()==(7381,h)
R=PolynomialRing(F,'q');ps=[R(list(v))for v in Q.rows()]
def coeff(p):v=p.list();return v+[0]*(prec-len(v))
A=Matrix(F,[coeff((ps[i]*ps[j]).truncate(prec))for i,j in pairs]);assert Tensor*Q2==A
B=Matrix(F,[coeff((ps[i-1]*ps[j-1]).truncate(prec))for i,j in T['basis_monomials']]);assert B==Q2
rows=(P/'ring12_tensor.txt').read_text().splitlines();assert rows[0]=='7381 360';assert Matrix(F,[list(map(int,s))for s in rows[1:]])==Tensor
print('C65_ALL_RING12_IDENTITIES_INDEPENDENT_PASS',round(time.monotonic()-t0,2),flush=True)
places=rd('rr_places_data.txt');cases=rd('rr_divisors_data.txt');answers=rd('rr12_results.txt');assert len(cases)==len(answers)==2

def conditions(ts,sign):
 selected=[(i,n)for i,n in ts if sign*n>0]
 return Matrix(F,[[c for i,n in selected for jet in places[i-1][2][j][:abs(n)]for c in jet]for j in range(g)])
def pidx(i,j):
 i,j=sorted([i,j]);return i*g-i*(i-1)//2+j-i
def mult(v):
 S=Matrix(F,g,len(pairs))
 for i,c in enumerate(v):
  if c:
   for j in range(g):S[j,pidx(i,j)]=1
 return S*Tensor
records=[]
for case,result in zip(cases,answers):
 mode,ii,ts=case;rm,ri,rts,degree,rp,rn,ff,hs,ranks,rankMF,rankfinal=result;assert[rm,ri,rts]==case
 positive=sum(places[i-1][0]*n for i,n in ts if n>0);negative=-sum(places[i-1][0]*n for i,n in ts if n<0);assert positive==negative==degree==36
 Cp=conditions(ts,1);Cm=conditions(ts,-1);assert Cp.ncols()==Cm.ncols()==degree and Cp.rank()==rp==Cm.rank()==rn==degree
 f=vector(F,ff);assert f and f*Cp==0;MF=mult(f);assert MF.rank()==rankMF==g;N=MF.right_kernel().basis_matrix().transpose();assert N.ncols()==h-g and MF*N==0;test=Cm
 for hh,oldrank in zip(hs,ranks):
  u=vector(F,hh);assert u*Cp==0;test=test.augment(mult(u)*N);assert test.rank()==oldrank
 assert test.rank()==rankfinal==g
 records.append(dict(source_degree=mode,variant=ii,necessary_point_variant=(ii==1),positive_degree=degree,negative_degree=degree,positive_vanishing_rank=degree,negative_vanishing_rank=degree,multiplication_by_F_rank=g,final_criterion_rank=g,H_sections=len(hs)))
 print('C65_CM6_NONPRINCIPAL_DIVISOR_INDEPENDENT_PASS',ii,degree,round(time.monotonic()-t0,2),flush=True)
out=dict(complete=True,level=65,prime=2,genus=g,canonical_q_precision=prec,weight4_Sturm_bound=672,ring_dimensions=[g,h],all_ordered_R1_products=g*g,all_symmetric_R1_products=len(pairs),all_q_product_identities_exact=True,source_basis_sha256=sha,nonprincipal_divisors=records,necessary_point_coefficients=[1,2,-3],extra_diagnostic_coefficients=[3,-2,-1])
if '--readonly'in sys.argv:assert json.loads((P/'ring_RR_sage_certificate.json').read_text())==out
else:(P/'ring_RR_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C65_CM6_RING_RR_INDEPENDENT_SAGE_COMPLETE_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
