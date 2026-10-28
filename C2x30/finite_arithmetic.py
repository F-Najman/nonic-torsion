"""Independent exact finite-field arithmetic for the eighteen source points.

Verifies primitive source levels and all54 actual quotient moduli maps.
Finite fields and elliptic addition are implemented here; Magma is not used.
The native RR/factorization stages are separately reproduced from scratch.
"""
import argparse,ast,re,runpy
from pathlib import Path
from sympy.polys.galoistools import gf_irreducible_p
from sympy.polys.domains import ZZ

def read_magma_data(path):
    s=path.read_text().replace('[*','[').replace('*]',']').replace(':=','=').replace('<','[').replace('>',']')
    s=re.sub(r'\btrue\b','True',s);s=re.sub(r'\bfalse\b','False',s)
    tree=ast.parse(s);out={}
    for stmt in tree.body:
        assert isinstance(stmt,ast.Assign) and len(stmt.targets)==1 and isinstance(stmt.targets[0],ast.Name)
        out[stmt.targets[0].id]=ast.literal_eval(stmt.value)
    return out

class Field:
    def __init__(self,pol):
        self.pol=tuple(x%7 for x in pol);self.n=len(pol)-1;assert self.pol[-1]==1
        assert gf_irreducible_p(list(reversed(self.pol)),7,ZZ)
        self.q=7**self.n;self.zero=Elt(self,[0]*self.n);self.one=Elt(self,[1]+[0]*(self.n-1));self.z=Elt(self,[0,1]+[0]*(self.n-2))
    def __call__(self,x):
        if isinstance(x,Elt):assert x.K is self;return x
        return Elt(self,[int(x)%7]+[0]*(self.n-1))

class Elt:
    def __init__(self,K,v):self.K=K;self.v=tuple(int(x)%7 for x in v);assert len(self.v)==K.n
    def __add__(self,b):
        b=self.K(b);return Elt(self.K,[(x+y)%7 for x,y in zip(self.v,b.v)])
    __radd__=__add__
    def __neg__(self):return Elt(self.K,[-x for x in self.v])
    def __sub__(self,b):return self+-self.K(b)
    def __rsub__(self,b):return self.K(b)+-self
    def __mul__(self,b):
        b=self.K(b);n=self.K.n;a=[0]*(2*n-1)
        for i,x in enumerate(self.v):
            if x:
                for j,y in enumerate(b.v):a[i+j]=(a[i+j]+x*y)%7
        for i in range(2*n-2,n-1,-1):
            c=a[i]
            if c:
                for j in range(n):a[i-n+j]=(a[i-n+j]-c*self.K.pol[j])%7
        return Elt(self.K,a[:n])
    __rmul__=__mul__
    def __pow__(self,n):
        if n<0:return self.inv()**(-n)
        r=self.K.one;b=self
        while n:
            if n&1:r=r*b
            b=b*b;n//=2
        return r
    def inv(self):assert self!=0;return self**(self.K.q-2)
    def __truediv__(self,b):return self*self.K(b).inv()
    def __rtruediv__(self,b):return self.K(b)*self.inv()
    def __eq__(self,b):
        if isinstance(b,int):return self.v==(b%7,)+(0,)*(self.K.n-1)
        return isinstance(b,Elt) and self.K is b.K and self.v==b.v

def evaluate(cs,x):
    out=x.K.zero
    for c in reversed(cs):out=out*x+c
    return out

class Curve:
    def __init__(self,K,b,c):self.K=K;self.a2=c;self.a4=b*(1-b-c);self.a6=K.zero
    def check(self,P):
        if P is None:return
        x,y=P;assert y*y==x*x*x+self.a2*x*x+self.a4*x+self.a6
    def add(self,P,Q):
        if P is None:return Q
        if Q is None:return P
        x,y=P;u,v=Q
        if x==u:
            if y+v==0:return None
            s=(3*x*x+2*self.a2*x+self.a4)/(2*y)
        else:s=(v-y)/(u-x)
        xx=s*s-self.a2-x-u;yy=-y+s*(x-xx);R=(xx,yy);self.check(R);return R
    def mul(self,n,P):
        assert n>=0;out=None
        while n:
            if n&1:out=self.add(out,P)
            P=self.add(P,P);n//=2
        return out
    def tate_xy(self,P):
        self.check(P);x,y=P;aa3=2*y;aa2=3*x+self.a2;aa4=3*x*x+2*x*self.a2+self.a4
        aaa1=2*aa4/aa3;aaa2=aa2-(aa4/aa3)**2
        b=-aaa2**3/aa3**2;c=-(aaa1*aaa2-aa3)/aa3
        r=b/c;s=c*c/(b-c);t=r*s-2*r+1
        return [(s-r)/t,t/(s*s-s-r+1)]

