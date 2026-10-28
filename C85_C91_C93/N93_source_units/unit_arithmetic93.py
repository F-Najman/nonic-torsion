"""Exact binary finite-field and short-series arithmetic for unit certificates."""
import ast,re

def lit(s,n):return ast.literal_eval(re.search(re.escape(n)+r':=(.*?);',s,re.S)[1].replace('[*','[').replace('*]',']').replace('<','(').replace('>',')'))
def matbits(s,n):
 m=re.search(re.escape(n)+r':=Matrix\(GF\(2\),(\d+),(\d+),(\[.*?\])\);',s,re.S);assert m,n;r,c=map(int,m.group(1,2));a=ast.literal_eval(m[3]);assert len(a)==r*c;return[sum(int(a[i*c+j])<<j for j in range(c))for i in range(r)]
def bits(a):return sum(int(v)<<i for i,v in enumerate(a))
def xorrows(v,rows):
 a=0
 while v:k=v&-v;a^=rows[k.bit_length()-1];v^=k
 return a
class Field:
 def __init__(self,p):
  self.p=p;self.d=p.bit_length()-1;self.squares=[self._rem(1<<(2*i))for i in range(self.d)]
 def _rem(self,a):
  while a.bit_length()>self.d:a^=self.p<<(a.bit_length()-1-self.d)
  return a
 def mul(self,a,b):
  c=0
  while b:
   if b&1:c^=a
   a<<=1;b>>=1
  return self._rem(c)
 def square(self,a):return xorrows(a,self.squares)
 def frob(self,a,n=1):
  for _ in range(n):a=self.square(a)
  return a
 def pow(self,a,n):
  b=1
  while n:
   if n&1:b=self.mul(b,a)
   a=self.square(a);n>>=1
  return b
 def inv(self,a):assert a;return self.pow(a,(1<<self.d)-2)
 def smul(self,a,b,n=None):
  n=len(a)if n is None else n;c=[0]*n
  for i,x in enumerate(a[:n]):
   if x:
    for j,y in enumerate(b[:n-i]):
     if y:c[i+j]^=self.mul(x,y)
  return c
 def sinv(self,a):
  n=len(a);iv=self.inv(a[0]);b=[iv]
  for k in range(1,n):
   z=0
   for i in range(1,k+1):z^=self.mul(a[i],b[k-i])
   b.append(self.mul(iv,z))
  return b
 def normalize(self,z):
  k=next(i for i,a in enumerate(z)if a[0]);iv=self.sinv(z[k]);return[self.smul(a,iv)for a in z]
 def key(self,z,d):
  zz=[a[0]if isinstance(a,list)else a for a in z];k=next(i for i,a in enumerate(zz)if a);iv=self.inv(zz[k]);zz=[self.mul(a,iv)for a in zz]
  divisors=[i for i in range(1,d)if d%i==0];gamma=next((a for a in zz if all(self.frob(a,i)!=a for i in divisors)),None)
  if gamma is None:gamma=next(a^b for a in zz for b in zz if all(self.frob(a^b,i)!=(a^b)for i in divisors))
  assert self.frob(gamma,d)==gamma;pp=[];v=1
  for _ in range(d):pp.append(v);v=self.mul(v,gamma)
  piv={}
  for i,a in enumerate(pp):
   c=1<<i
   while a:
    k=a.bit_length()-1
    if k in piv:b,e=piv[k];a^=b;c^=e
    else:piv[k]=(a,c);break
   assert a
  def coords(a):
   c=0
   while a:
    k=a.bit_length()-1;assert k in piv;b,e=piv[k];a^=b;c^=e
   return c
  pol=(1<<d)|coords(v);return pol,tuple(coords(a)for a in zz)

def sadd(*rows):
 c=[0]*len(rows[0])
 for a in rows:c=[x^y for x,y in zip(c,a)]
 return c

def unit_series(F,a2,a6,xx,yy):
 n=len(xx);one=[1]+[0]*(n-1);aa2=sadd(xx,[a2]+[0]*(n-1));w=F.smul(sadd(F.smul(xx,xx),yy),F.sinv(xx));aa=sadd(aa2,F.smul(w,w),w);b=F.smul(F.smul(F.smul(aa,aa),aa),F.sinv(F.smul(xx,xx)));c=F.smul(sadd(aa,xx),F.sinv(xx));r=F.smul(b,F.sinv(c));ss=F.smul(F.smul(c,c),F.sinv(sadd(b,c)));t=sadd(F.smul(r,ss),one);return F.smul(sadd(ss,r),F.sinv(t)),F.smul(t,F.sinv(sadd(F.smul(ss,ss),ss,r,one)))

def psi93(F,xx,a6):
 n=len(xx);one=[1]+[0]*(n-1);zero=[0]*n;xx2=F.smul(xx,xx);xx3=F.smul(xx2,xx);xx4=F.smul(xx2,xx2);xx6=F.smul(xx3,xx3);ps=[zero,one,xx,sadd(xx4,xx3,a6),sadd(xx6,F.smul(a6,xx2))];iv=F.sinv(xx)
 for k in range(5,94):
  m=k//2
  if k%2:v=sadd(F.smul(ps[m+2],F.smul(F.smul(ps[m],ps[m]),ps[m])),F.smul(ps[m-1],F.smul(F.smul(ps[m+1],ps[m+1]),ps[m+1])))
  else:v=F.smul(F.smul(ps[m],iv),sadd(F.smul(ps[m+2],F.smul(ps[m-1],ps[m-1])),F.smul(ps[m-2],F.smul(ps[m+1],ps[m+1]))))
  ps.append(v)
 return ps[93]

def decode_record(rec):
 p,aa,pp,zz,xy,jt=rec;F=Field(bits(p));return{'F':F,'field_poly':p,'ainv':[bits(a)for a in aa],'point':[bits(a)for a in pp],'z':[bits(a)for a in zz],'units':[bits(a)for a in xy],'xjet':[bits(a)for a in jt[0]],'yjet':[bits(a)for a in jt[1]],'jets':[[bits(a)for a in row]for row in jt[2]]}

def verify_record(rec,RX,RY):
 v=decode_record(rec);F=v['F'];aa=v['ainv'];xx,yy=v['xjet'],v['yjet'];zz=v['jets'];assert F.d in[9,18]and len(zz)==261 and all(len(a)==5 for a in zz);assert aa[0]==1 and aa[2]==aa[3]==0 and xx[0]==v['point'][0]and yy[0]==v['point'][1]and[a[0]for a in zz]==v['z'];assert any(a[1]for a in zz)
 a6=[aa[4],1,0,0,0];assert not any(psi93(F,xx,a6));de=psi93(F,[xx[0],1,0,0,0],[aa[4],0,0,0,0])[1];assert de
 curve=sadd(F.smul(yy,yy),F.smul(xx,yy),F.smul(F.smul(xx,xx),xx),[F.mul(aa[1],a)for a in F.smul(xx,xx)],a6);assert not any(curve)
 uv=unit_series(F,aa[1],a6,xx,yy);assert[a[0]for a in uv]==v['units']
 for rels,u in zip([RX,RY],uv):
  for row in rels:
   b=[xorrows(row&((1<<261)-1),[a[n]for a in zz])for n in range(5)];a=[xorrows(row>>261,[a[n]for a in zz])for n in range(5)];assert not any(sadd(b,F.smul(u,a)))
 v['unitjets']=uv;return v

def transform(F,z,D):
 cols=[[a[n]for a in z]for n in range(len(z[0]))];return F.normalize([[xorrows(row,col)for col in cols]for row in D])
