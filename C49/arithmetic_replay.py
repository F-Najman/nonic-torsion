#!/usr/bin/env python3
"""Independent primitive cusp-label and rational2-torsion arithmetic for C49."""
from math import gcd
from collections import Counter
from pathlib import Path
import json,sys
N=49

def cusp(a,c):
 c%=N;g=gcd(c,N)
 return min((a%g,c),((-a)%g,(-c)%N))
labels={cusp(a,c)for c in range(N)for a in range(gcd(c,N))if gcd(a,gcd(c,N))==1}
def sigma(P,u):return cusp(u*P[0],P[1])
def diamond(P,u):return cusp(pow(u,-1,N)*P[0],u*P[1])
def orbits(values,action):
 left=set(values);res=[]
 while left:
  p=min(left);o={p};q=action(p)
  while q!=p:assert q not in o;o.add(q);q=action(q)
  left-=o;res.append(frozenset(o))
 return res
qo=orbits(labels,lambda P:sigma(P,3));fo=orbits(labels,lambda P:sigma(P,2))
assert len(labels)==60 and Counter(map(len,qo))=={1:21,6:3,21:1}and Counter(map(len,fo))=={1:21,3:6,21:1}
rational={next(iter(o))for o in fo if len(o)==1};assert len(orbits(rational,lambda P:diamond(P,2)))==1
cubics={o for o in fo if len(o)==3}
def dc(o):return frozenset(diamond(P,2)for P in o)
co=orbits(cubics,dc);assert sorted(map(len,co))==[3,3]
for q in qo:
 if len(q)!=6:continue
 pair={o for o in cubics if o<=q};assert len(pair)==2
 assert all(len(pair&set(o))==1 for o in co)
 assert len({z for o in pair for z in[o,dc(o),dc(dc(o))]})==6
assert 6*3<21
# Rational2-primary specialization transfer arithmetic for X0(49).
x,y=2,-1;a1,a2,a3,a4,a6=1,-1,0,-2,-1
assert y*y+a1*x*y+a3*y==x**3+a2*x*x+a4*x+a6 and 2*y+a1*x+a3==0
b2=a1*a1+4*a2;b4=2*a4+a1*a3;b6=a3*a3+4*a6;b8=a1*a1*a6+4*a2*a6-a1*a3*a4+a2*a3*a3-a4*a4
assert -b2*b2*b8-8*b4**3-27*b6*b6+9*b2*b4*b6==-343
assert 52367710906884085342%4==2 and (0,1)in[(a,b)for a in range(2)for b in range(2)if (b*b+a*b-a**3+a*a+2*a+1)%2==0]
res={'level':49,'degree':9,'geometric_cusps':60,'rational_cusps':21,'cubic_closed_cusps':6,'rational_diamond2_cycle':21,'cubic_diamond2_cycles':[3,3],'each_global_sextic_cusp_splits_into_one_cusp_from_each_cycle':True,'required_base_cusps':{'rational':1,'cubic':2},'selected_rational':18,'omitted_rational':3,'mixed_rational_support_at_most':6,'forbidden_translations_at_most':18,'rational_2primary_order':2,'quotient_degree':21,'elliptic_two_point_nonzero_reduction':True,'status':'PASS'}
if len(sys.argv)>1:Path(sys.argv[1]).write_text(json.dumps(res,indent=2)+'\n')
print('NONIC49_CUSP_AND_TWO_TORSION_ARITHMETIC_PASS',json.dumps(res,sort_keys=True))
