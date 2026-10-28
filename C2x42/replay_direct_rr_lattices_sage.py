from pathlib import Path
import ast,json,sys,time,resource,os
os.nice(10);resource.setrlimit(resource.RLIMIT_CPU,(1770,1800));resource.setrlimit(resource.RLIMIT_AS,(10*10**9,10*10**9))
from sage.all import GF,PolynomialRing,matrix,identity_matrix,lcm,PowerSeriesRing,divisors
P=Path(__file__).resolve().parent/'data';t0=time.monotonic()
def rd(p):return ast.literal_eval(p.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False'))
models=[rd(P/d/'model_bases.txt')for d in['cyclic42_A11_seven','quartic_alt_A11','CM_cyclic42_lattice']];assert models[0]==models[1]==models[2]
n,fp,coords,bf,bi=models[0];assert n==15
k=GF(5);R=PolynomialRing(k,'x');K=R.fraction_field();x=K.gen();S=PolynomialRing(K,'y');y=S.gen()
def cfun(c):return K(R(c[0]))/K(R(c[1]))
def mat(b):return matrix(K,[[cfun(c)for c in row]for row in b])
f=S([cfun(c)for c in fp]).monic();assert f.degree()==n
assert [cfun(c)for c in coords[0]]==[x]+[K.zero()]*(n-1)and [cfun(c)for c in coords[1]]==[K.zero(),K.one()]+[K.zero()]*(n-2)
BF=mat(bf);BI=mat(bi);BFI=BF.inverse();BII=BI.inverse();SP=PolynomialRing(k,2,names=('X','Y'));X,Y=SP.gens();assert all(c[1]==[1]for c in fp);Fraw=sum(sum(SP(v)*X**j for j,v in enumerate(c[0]))*Y**i for i,c in enumerate(fp));res=Fraw.resultant(Fraw.derivative(Y),Y);rr=R(res.univariate_polynomial());lc=R(fp[-1][0]);disc=(-1)**(n*(n-1)//2)*K(rr)/K(lc)**(2*n-1);assert disc!=0;print('C242_GENUS25_FAST_DISCRIMINANT_READY',time.monotonic()-t0,flush=True);df=disc*BF.det()**2;di=disc*BI.det()**2;assert df.denominator().degree()==0
fin=int(df.numerator().degree());inf=int(di.denominator().degree()-di.numerator().degree());assert fin>=0 and inf>=0 and fin+inf==2*25-2+2*n==78
for name,B in [('finite',BF),('infinite',BI)]:
 for j,row in enumerate(B.rows()):
  a=S(list(row));T=matrix(K,[[(a*y**i%f)[h]for h in range(n)]for i in range(n)]);cp=T.charpoly()
  if name=='finite':assert all(c.denominator().degree()==0 for c in cp)
  else:assert all(not c or c.numerator().degree()<=c.denominator().degree()for c in cp)
  print('C242_GENUS25_BASIS_INTEGRAL',name,j+1,round(time.monotonic()-t0,2),flush=True)
orderdata=dict(complete=True,degree=15,genus=25,finite_discriminant_degree=fin,infinite_discriminant_valuation=inf,total_different_degree=78,independent_integrality_characteristic_polynomials=30)
if '--readonly'not in sys.argv:(P/'direct_rr_order_certificate.json').write_text(json.dumps(orderdata,indent=2)+'\n')
else:assert json.loads((P/'direct_rr_order_certificate.json').read_text())==orderdata
geo=json.loads((P/'final_divisor_geometry.json').read_text());assert geo['complete'];points=geo['points'];cases=geo['cases'];cache={};R0=PolynomialRing(k,'z')
def pt(key):
 if key not in cache:
  d=points[key];L=GF(5**d['degree'],'z',modulus=R0(d['polynomial']));z=L.gen();xx=sum(L(a)*z**i for i,a in enumerate(d['x_coefficients']));yy=sum(L(a)*z**i for i,a in enumerate(d['y_coefficients']));base=tuple(map(int,xx.minpoly()));dx=len(base)-1;assert d['degree']%dx==0
  cache[key]=(L,xx,yy,base,dx)
 return cache[key]
def ce(c,z):
 return sum(z.parent()(a)*z**i for i,a in enumerate(c.numerator()))/sum(z.parent()(a)*z**i for i,a in enumerate(c.denominator()))
def pol_in_series(a,zs,RY):return RY([ce(c,zs)for c in a])
def series_root(xx,yy,prec):
 L=xx.parent();PS=PowerSeriesRing(L,'e',default_prec=prec);e=PS.gen();xs=(PS(xx)+e).add_bigoh(prec);cs=[ce(c,xs)for c in f];ys=PS(yy).add_bigoh(prec)
 for _ in range(7):
  Fv=sum(c*ys**j for j,c in enumerate(cs));dv=sum(j*c*ys**(j-1)for j,c in enumerate(cs)if j);assert dv[0]!=0;ys=(ys-Fv/dv).add_bigoh(prec)
 assert sum(c*ys**j for j,c in enumerate(cs)).valuation()>=prec
 return xs,ys
# Every nonprimitive-x point is quadratic over its base x-field; all other points use x itself.
assert sum(not p['x_primitive']for p in points.values())==1

def validate_ideal(II,D,label):
 pos={q:m for q,m in D.items()if m>0};neg={q:-m for q,m in D.items()if m<0};groups={}
 for q,m in pos.items():groups.setdefault(pt(q)[3],[]).append((q,m))
 Tp=R.one()
 for bp,qs in groups.items():Tp*=R(bp)**max(m for _,m in qs)
 C=Tp*II*BFI;assert all(a.denominator().degree()==0 for a in C.list());dd=C.det();assert dd.denominator().degree()==0 and dd.numerator().degree()==n*Tp.degree()
 # Local equation orders coincide with the already certified maximal order at these fibres.
 for q in D:
  L,xx,yy,bp,dx=pt(q);assert f.leading_coefficient()==1
  assert all(sum(L(a)*xx**i for i,a in enumerate(c.denominator()))!=0 for c in f)
  det=BF.det();assert ce(det,xx)!=0 and sum(L(a)*xx**i for i,a in enumerate(det.denominator()))!=0
  FY=sum(j*ce(c,xx)*yy**(j-1)for j,c in enumerate(f)if j);assert FY!=0
 polTI=[S([Tp*c for c in row])for row in II.rows()]
 # At a positive base fibre, retain precisely the allowed branches, with all other branches killed.
 for bp,qs in groups.items():
  qmax=max(qs,key=lambda q:points[q[0]]['degree'])[0];L,xx,yy,_,dx=pt(qmax);RP=PolynomialRing(L,'Y');Y=RP.gen();selected=RP.one();targets=[]
  for q,m in qs:
   LL,x0,y0,_,dxx=pt(q);assert dxx==dx
   if LL is L:assert x0==xx;yval=y0
   else:
    assert points[q]['x_primitive'];em=LL.hom([xx],L);assert em(x0)==xx;yval=em(y0)
   h=points[q]['degree']//dx;fac=RP.one()
   for j in range(h):fac*=Y-yval**(5**(dx*j))
   selected*=fac;targets.append((yval,h,m))
  ff=RP([ce(c,xx)for c in f]);assert ff%selected==0 and selected.is_squarefree();other=ff.quo_rem(selected)[0];assert selected.gcd(other)==1
  order=max(m for _,m in qs);PS=PowerSeriesRing(L,'e',default_prec=order);ep=PS.gen();xs=(PS(xx)+ep).add_bigoh(order);RY=PolynomialRing(PS,'Y');Fs=pol_in_series(f,xs,RY);As=RY([PS(a)for a in selected]);Bs=RY([PS(a)for a in other]);_,inv,_=selected.xgcd(other)
  for r in range(1,order):
   err=Fs-As*Bs;Er=RP([c[r]for c in err]);db=(Er*inv)%other;da=(Er-selected*db).quo_rem(other);assert da[1]==0;As+=ep**r*RY([PS(a)for a in da[0]]);Bs+=ep**r*RY([PS(a)for a in db])
  assert all(c.valuation()>=order for c in Fs-As*Bs)
  for a in polTI:
   gs=pol_in_series(a,xs,RY);rem=gs.quo_rem(Bs)[1];assert all(c.valuation()>=order for c in rem)
  for yy,h,m in targets:
   required=order-m
   if not required:continue
   for j in range(h):
    xs,ys=series_root(xx,yy**(5**(dx*j)),required)
    for a in polTI:assert sum(ce(c,xs)*ys**i for i,c in enumerate(a)).valuation()>=required
 # At each negative point, exact ordinary jets impose the additional vanishing order.
 for q,m in neg.items():
  L,xx,yy,bp,dx=pt(q);order=max([v for _,v in groups.get(bp,[])],default=0)+m;xs,ys=series_root(xx,yy,order)
  for a in polTI:assert sum(ce(c,xs)*ys**i for i,c in enumerate(a)).valuation()>=order
 assert sum(points[q]['degree']*m for q,m in D.items())==0
 print('C242_GENUS25_EXACT_FINITE_IDEAL_CERTIFIED',label,'norm',n*Tp.degree(),round(time.monotonic()-t0,2),flush=True)
 return int(n*Tp.degree())

def popov(II,JJ,expected,label):
 M=II*JJ.inverse();den=lcm([a.denominator()for a in M.list()]);Mat=matrix(R,n,[(den*a).numerator()for a in M.list()]);orig=matrix(R,n,n,list(Mat.list()));U=identity_matrix(R,n)
 def piv(row):return max((int(a.degree()),j)for j,a in enumerate(row)if a)
 steps=0
 while True:
  ps=[piv(row)for row in Mat.rows()];pair=next(((i,j)for i in range(n)for j in range(i)if ps[i][1]==ps[j][1]),None)
  if pair is None:break
  i,j=pair
  if ps[i][0]<ps[j][0]:i,j=j,i
  c=ps[i][1];s=-Mat[i,c].leading_coefficient()/Mat[j,c].leading_coefficient()*R.gen()**(ps[i][0]-ps[j][0]);Mat.add_multiple_of_row(i,j,s);U.add_multiple_of_row(i,j,s);steps+=1
 # Each recorded update is a polynomial elementary row addition, hence U is unimodular.
 assert Mat==U*orig
 ds=[piv(r)[0]for r in Mat.rows()];dim=sum(max(0,den.degree()-d+1)for d in ds);assert dim==expected
 print('C242_GENUS25_INDEPENDENT_RR_PASS',label,dim,round(time.monotonic()-t0,2),flush=True)
 return dict(dimension=int(dim),denominator_degree=int(den.degree()),row_degrees=ds,unimodular_row_additions=steps)

answers={}
for label,D in cases.items():
 if label.startswith('seven_'):
  i=int(label[6:]);data=rd(P/'cyclic42_A11_seven'/f'ideal_bases_case{i}.txt');I,J=data[1:];expected=1 if i==1 else 0
 elif label.startswith('alt_'):
  i=int(label[4:]);I,J=rd(P/'quartic_alt_A11'/f'ideal_bases_map{i}.txt');expected=0
 else:
  i=int(label[3:]);data=rd(P/'CM_cyclic42_lattice'/f'ideal_bases_case{i}.txt');assert data[0]==i;I,J=data[1:];expected=0
 assert J==bi;II=mat(I);JJ=mat(J);norm=validate_ideal(II,D,label);answer=popov(II,JJ,expected,label);answer['finite_ideal_norm_index']=norm;answers[label]=answer
out=dict(complete=True,prime=5,genus=25,field_degree=15,integral_bases=orderdata,all_14_exact_RR_ideals_verified=True,all_14_RR_dimensions_independently_verified=True,answers=answers)
if '--readonly'in sys.argv:assert json.loads((P/'direct_rr_lattice_certificate.json').read_text())==out
else:(P/'direct_rr_lattice_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
print('C242_ALL_FINAL_RR_LATTICES_INDEPENDENT_SAGE_PASS',json.dumps(out),'seconds',time.monotonic()-t0,flush=True)
