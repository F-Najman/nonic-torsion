from pathlib import Path
import ast,re,json,time,sys
from sage.all import GF,PolynomialRing,LaurentSeriesRing,Matrix
P=Path(__file__).resolve().parents[1]/'data'/'source85';t0=time.monotonic()
def lit(s):return ast.literal_eval(s.replace('[*','[').replace('*]',']'))
def mat(s,name):
 m=re.search(name+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,name
 return Matrix(GF(2),int(m[1]),int(m[2]),lit(m[3]))
raw=lit((P/'quadratic_cusp_arc_data.txt').read_text());width,pol,oldprec,pivot,inds,idx,zco,xco,yco=raw;assert width==3 and pivot==idx==1
s=(P/'canonical_unit_data.m').read_text();RX=mat(s,'UnitXRelations');RY=mat(s,'UnitYRelations');assert RX.nrows()==49 and RY.nrows()==47
F=GF(4,'a',modulus=PolynomialRing(GF(2),'x')(pol));a=F.gen();prec=220;LS=LaurentSeriesRing(F,'q',default_prec=prec);q=LS.gen();w=19;u=a*q;s1=LS(0);a6=LS(0);xp=u/(1+u)**2;yp=u*u/(1+u)**3
for n in range(1,prec//w+2):
 if w*n+1<prec:v=a*q**(w*n+1);xp+=v/(1+v)**2;yp+=v*v/(1+v)**3
 if w*n-1<prec:v=a**-1*q**(w*n-1);xp+=v/(1+v)**2;yp+=v/(1+v)**3
 if w*n<prec:v=q**(w*n);s1+=F(n)*v/(1+v);a6+=F((5*n**3+7*n**5)//12)*v/(1+v)
yp+=s1;assert(yp*yp+xp*yp+xp**3+s1*xp+a6).valuation()>=prec
A4=xp*xp+s1+yp;A2=xp+(A4/xp)**2+A4/xp;b=A2**3/xp**2;c=(A2+xp)/xp;r=b/c;ss=c*c/(b+c);tt=r*ss+1;xx=(ss+r)/tt;yy=tt/(ss*ss+ss+r+1)
def decode(v):return sum(F(c)*a**j for j,c in enumerate(v))
assert all(xx[n]==decode(v)for n,v in zip(range(-1,31),xco))and all(yy[n]==decode(v)for n,v in zip(range(-1,31),yco))
M=[[LS(int(v[j]))+zz*int(v[85+j])for j in range(85)]for R,zz in[(RX,xx),(RY,yy)]for v in R.rows()];cols=[j for j in range(85)if j!=idx-1]
A=[[q*M[i-1][j]for j in cols]for i in inds];rhs=[-q*M[i-1][idx-1]for i in inds];n=84;perm=list(range(n));detval=0
for k in range(n):
 val,i,j=min((A[i][j].valuation(),i,j)for i in range(k,n)for j in range(k,n));assert val<A[i][j].precision_absolute();detval+=int(val)
 A[k],A[i]=A[i],A[k];rhs[k],rhs[i]=rhs[i],rhs[k]
 if j!=k:
  for row in A:row[k],row[j]=row[j],row[k]
  perm[k],perm[j]=perm[j],perm[k]
 piv=A[k][k]
 for j in range(k+1,n):A[k][j]/=piv
 rhs[k]/=piv;A[k][k]=LS(1)
 for i in range(k+1,n):
  v=A[i][k]
  if v:
   for j in range(k+1,n):A[i][j]-=v*A[k][j]
   rhs[i]-=v*rhs[k]
  A[i][k]=LS(0)
 if(k+1)%20==0:print('C57_QUADRATIC_ARC_INDEPENDENT_PIVOT',k+1,detval,round(time.monotonic()-t0,2),flush=True)
assert detval==95
y=[LS(0)]*n
for k in reversed(range(n)):y[k]=rhs[k]-sum(A[k][j]*y[j]for j in range(k+1,n))
sol=[LS(0)]*n
for k in range(n):sol[perm[k]]=y[k]
z=[LS(1)]+sol;assert all(v.valuation()>=0 and v.precision_absolute()>=20 for v in z)
assert all(z[i][j]==decode(zco[i][j])for i in range(85)for j in range(20));assert all(sum(v*w for v,w in zip(row,z)).valuation()>=20 for row in M);assert any(v[1]for v in z)
out=dict(complete=True,level=57,prime=2,genus=85,source_Tate_width=19,dual_j_width=3,toric_root_order=3,independent_Tate_precision=prec,canonical_coordinates=85,all_unit_relations=96,normalized_coordinate=1,independent_selected_minor_valuation=95,verified_canonical_jets=20,nonzero_tangent=True,genuine_Tate_family=True)
if '--readonly'in sys.argv:assert json.loads((P/'quadratic_arc_sage_certificate.json').read_text())==out
else:(P/'quadratic_arc_sage_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C57_QUADRATIC_CUSP_ARC_INDEPENDENT_SAGE_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
