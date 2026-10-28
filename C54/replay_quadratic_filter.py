#!/usr/bin/env python3
"""Independent prime-polynomial arithmetic for every quadratic skip certificate."""
from pathlib import Path
import ast,json,sys
P=Path(__file__).resolve().parent/'data'

def main():
 D=json.loads((P/'quadratic_prepass_results.json').read_text());F=json.loads((P/'quadratic_prepass_function.json').read_text());assert D['complete'] and F['complete']
 assert F['function_divisor_even'] and F['both_orders_exact'] and F['H5_odd_order']==472473
 raw=ast.literal_eval((P/'good9_pointdata.txt').read_text().replace('[*','[').replace('*]',']'));pol,records=raw
 assert pol==D['field_polynomial'] and len(pol)==10 and pol[-1]%5==1 and len(records)==16135
 n=9;q=5**n;zero=(0,)*n;one=(1,)+(0,)*(n-1)
 def elt(a):return tuple([x%5 for x in a]+[0]*(n-len(a)))
 def add(a,b):return tuple((x+y)%5 for x,y in zip(a,b))
 def mul(a,b):
  c=[0]*17
  for i,x in enumerate(a):
   if x:
    for j,y in enumerate(b):c[i+j]+=x*y
  for k in range(16,8,-1):
   t=c[k]%5
   if t:
    for j in range(9):c[k-9+j]-=t*pol[j]
  return tuple(c[i]%5 for i in range(9))
 def power(a,e):
  out=one
  while e:
   if e&1:out=mul(out,a)
   a=mul(a,a);e>>=1
  return out
 def evalpoly(cs,a):
  out=zero
  for c in reversed(cs):out=add(mul(out,a),elt([c]))
  return out
 fc=int(F['base_cusp_value']);assert fc in range(1,5);invfc=pow(fc,-1,5)
 coeff=F['filter_function'];signs=[]
 for i,r in enumerate(records,1):
  x,y=map(elt,r[0]);value=zero;ypow=one;inverses={};fallback=False
  for num,den in coeff:
   de=evalpoly(den,x)
   if de==zero:fallback=True;break
   if de not in inverses:inverses[de]=power(de,q-2)
   value=add(value,mul(mul(evalpoly(num,x),inverses[de]),ypow));ypow=mul(ypow,y)
  if fallback or value==zero:sign=0
  else:
   value=tuple(a*invfc%5 for a in value);v=power(value,(q-1)//2)
   assert v in (one,elt([4]));sign=1 if v==one else -1
  signs.append(sign)
 assert signs==D['signs']
 skips=[i for i,v in enumerate(signs,1) if v==-1];remaining=[i for i,v in enumerate(signs,1) if v!=-1];fallback=[i for i,v in enumerate(signs,1) if v==0]
 assert skips==D['nonsquare_indices'] and remaining==D['remaining_indices'] and fallback==D['fallback_indices']
 print('NONIC54_QUADRATIC_INTEGER_REPLAY '+json.dumps({'records':len(records),'nonsquare_certificates':len(skips),'remaining':len(remaining),'fallback':len(fallback),'all_finite_field_arithmetic_independently_recomputed':True,'complete':True}))
if __name__=='__main__':main()
