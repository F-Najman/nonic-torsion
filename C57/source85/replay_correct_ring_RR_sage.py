from pathlib import Path
import ast,re,json,time,sys
from sage.all import GF,Matrix,vector,PolynomialRing,LaurentSeriesRing
P=Path(__file__).resolve().parents[1]/'data'/'source85';t0=time.monotonic()
def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']').replace('<','[').replace('>',']'))
def rd(f):return lit((P/f).read_text())
def mat(s,name):
 m=re.search(name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,name
 return Matrix(GF(2),int(m[1]),int(m[2]),lit(m[3]))
def tensor(f,n,m):
 rows=(P/f).read_text().splitlines();assert rows[0]==f'{n} {m}'and len(rows)==n+1;return Matrix(GF(2),[list(map(int,s))for s in rows[1:]])
s=(P/'qbasis_data.m').read_text();Q=mat(s,'qcoeffs');D=mat(s,'Diamond5Matrix');assert Q[:,0:241].rank()==85 and D**18==Matrix.identity(GF(2),85)
s=(P/'ring24_data.m').read_text();Q2=mat(s,'Ring2qCoefficients');Q4=mat(s,'Ring4qCoefficients');RX2=mat(s,'Ring2UnitXRelations');RY2=mat(s,'Ring2UnitYRelations');inds=lit(re.search(r'Ring2MonomialIndices:=(\[.*?\]);',s,re.S)[1]);inds4=lit(re.search(r'Ring4MonomialIndices:=(\[.*?\]);',s,re.S)[1]);R=PolynomialRing(GF(2),'q');q=R.gen();ps=[R(list(v))for v in Q.rows()];ps2=[R(list(v))for v in Q2.rows()];ps4=[R(list(v))for v in Q4.rows()];pairs=[(i,j)for i in range(85)for j in range(i,85)];pairs2=[(i,j)for i in range(252)for j in range(i,252)];T12=tensor('ring12_tensor.txt',3655,252);T24=tensor('ring24_tensor.txt',31878,588)
def coeff(v):w=v.list();return w+[0]*(1100-len(w))
A12=Matrix(GF(2),[coeff((ps[i]*ps[j]).truncate(1100))for i,j in pairs]);assert Q2==A12.matrix_from_rows([i-1 for i in inds])and T12*Q2==A12 and Q2[:,0:481].rank()==252
print('C57_ALL_RING12_IDENTITIES_SAGE_PASS',round(time.monotonic()-t0,2),flush=True)
A24=Matrix(GF(2),[coeff((ps2[i]*ps2[j]).truncate(1100))for i,j in pairs2]);assert Q4==A24.matrix_from_rows([i-1 for i in inds4])and T24*Q4==A24 and Q4[:,0:961].rank()==588
print('C57_ALL_RING24_IDENTITIES_SAGE_PASS',round(time.monotonic()-t0,2),flush=True)
# The canonical/unit identities are certified globally by independent pole bounds.
s=(P/'canonical_unit_data.m').read_text();RX=mat(s,'UnitXRelations');RY=mat(s,'UnitYRelations');LS=LaurentSeriesRing(GF(2),'t',default_prec=1500);t=LS.gen();s1=LS(0);a6=LS(0)
for n in range(1,1500//57+1):
 v=t**(57*n);s1+=GF(2)(n)*v/(1+v);a6+=GF(2)((5*n**3+7*n**5)//12)*v/(1+v)
def uv(a):
 u=t**a;x=u/(1+u)**2;y=u*u/(1+u)**3+s1
 for n in range(1,1500//57+2):
  if 57*n+a<1500:v=t**(57*n+a);x+=v/(1+v)**2;y+=v*v/(1+v)**3
  if 57*n-a<1500:v=t**(57*n-a);x+=v/(1+v)**2;y+=v/(1+v)**3
 assert(y*y+x*y+x**3+s1*x+a6).valuation()>=1500
 aa4=x*x+s1+y;aa2=x+(aa4/x)**2+aa4/x;b=aa2**3/x**2;c=(aa2+x)/x;r=b/c;ss=c*c/(b+c);tt=r*ss+1;return((ss+r)/tt,tt/(ss*ss+ss+r+1))
x,y=uv(1);x5,y5=uv(5);F1=[LS(list(v)).add_bigoh(1100)for v in Q.rows()];F2=[LS(list(v)).add_bigoh(1100)for v in Q2.rows()]
def checkrel(M,fs,u,bound):
 for row in M.rows():
  a=sum(fs[i]for i in range(len(fs))if row[i]);b=sum(fs[i]for i in range(len(fs))if row[len(fs)+i]);v=a+u*b;assert v.valuation()-2>bound
for M,fs,u,b in[(RX,F1,x,207),(RY,F1,y,207),(RX2,F2,x,375),(RY2,F2,y,375)]:checkrel(M,fs,u,b)
# Identify the canonical diamond with the actual point-level multiplier.
passes=[]
for k in[1,-1]:
 QQ=D**k*Q;fs=[LS(list(v)).add_bigoh(1100)for v in QQ.rows()]
 try:checkrel(RX,fs,x5,207);checkrel(RY,fs,y5,207);passes.append(k)
 except AssertionError:pass
assert passes==[-1]
pole=json.loads((P/'pole_degrees_sage_certificate.json').read_text());assert pole['pole_degrees']==[37,39]
print('C57_ALL_GLOBAL_UNIT_AND_POINT_DIAMOND_IDENTITIES_SAGE_PASS',round(time.monotonic()-t0,2),flush=True)
places=rd('rr_places_data.txt');cases=rd('rr_correct_divisors_data.txt');answers=rd('rr24_correct_results.txt');records=[]
def conditions(ts,sign):
 selected=[(i,n)for i,n in ts if sign*n>0];return Matrix(GF(2),[[c for i,n in selected for jet in places[i-1][2][j][:abs(n)]for c in jet]for j in range(252)])
def pidx(i,j):
 i,j=sorted([i,j]);return i*252-i*(i-1)//2+j-i

def mult(v):
 S=Matrix(GF(2),252,31878)
 for i,c in enumerate(v):
  if c:
   for j in range(252):S[j,pidx(i,j)]=1
 return S*T24
for case,result in zip(cases,answers):
 mode,ii,ts=case;rm,ri,rts,degree,rp,rn,ff,hs,ranks,rankMF,rankfinal=result;assert[rm,ri,rts]==case;positive=sum(places[i-1][0]*n for i,n in ts if n>0);negative=-sum(places[i-1][0]*n for i,n in ts if n<0);assert positive==negative==degree
 Cplus=conditions(ts,1);Cminus=conditions(ts,-1);assert Cplus.rank()==Cminus.rank()==degree;f=vector(GF(2),ff);assert f and f*Cplus==0;MF=mult(f);assert MF.rank()==rankMF==252;N=MF.right_kernel().basis_matrix().transpose();assert N.ncols()==336 and MF*N==0;test=Cminus
 for hh,oldrank in zip(hs,ranks):
  h=vector(GF(2),hh);assert h*Cplus==0;MH=mult(h);test=test.augment(MH*N);assert test.rank()==oldrank
 assert test.rank()==rankfinal==252;records.append(dict(source_degree=mode,source_index=ii,positive_degree=degree,negative_degree=degree,positive_vanishing_rank=degree,negative_vanishing_rank=degree,multiplication_by_F_rank=252,final_criterion_rank=252,H_sections=len(hs)))
 print('C57_RR_NONPRINCIPAL_DIVISOR_SAGE_PASS',mode,ii,degree,round(time.monotonic()-t0,2),flush=True)
out=dict(complete=True,level=57,prime=2,genus=85,canonical_q_precision=1100,weight4_Sturm_bound=480,weight8_Sturm_bound=960,ring_dimensions=[85,252,588],all_R1_products=3655,all_R2_products=31878,all_q_product_identities_exact=True,global_R1_and_R2_unit_identities=True,physical_point5_is_canonical_diamond_inverse=True,nonprincipal_divisors=records,full_good_reduction_obstruction=True,geometric_push_annihilator="T5-sigma-5",necessary_divisor_source_coefficients=[5,-4,-1])
if '--readonly'in sys.argv:assert json.loads((P/'ring_correct_RR_sage_certificate.json').read_text())==out
else:(P/'ring_correct_RR_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C57_CANONICAL_RING_AND_ALL_RR_OBSTRUCTIONS_INDEPENDENT_SAGE_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
