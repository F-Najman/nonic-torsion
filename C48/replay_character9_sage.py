#!/usr/bin/env python3
"""Independent finite-field replay using Sage, without Magma or class maps."""
import os,resource
os.nice(10)
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
resource.setrlimit(resource.RLIMIT_CPU,(600,600))
resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))
from pathlib import Path
import ast,json,time
from sage.all import GF,PolynomialRing
P=Path(__file__).resolve().parent

def literal(f):
 s=f.read_text().replace('[*','[').replace('*]',']')
 return {n.targets[0].id:ast.literal_eval(n.value) for n in ast.parse(s).body}
s=(P/'data/parameter9/records.m').read_text();poly=ast.literal_eval(s.split('FieldPolynomial:=',1)[1].split(';',1)[0]);records=ast.literal_eval(s.split('Records:=',1)[1].split(';\nComplete',1)[0].replace('[*','[').replace('*]',']'))
D=literal(P/'data/finite_characters_v2/function_data.py');F=literal(P/'data/character_filter9/results.py');geo=literal(P/'data/degree8/clean_run/factorbase_geometry.py')
k=GF(5);rk=PolynomialRing(k,'w');K=GF(5**9,'z',modulus=rk(poly));z=K.gen();R=PolynomialRing(K,'v');one=K.one();zero=K.zero()
def elt(raw):return K(rk(raw))
def rationaldata(fn):return [(R(v[0]),R(v[1])) for v in fn]
functions={i:rationaldata(D['functions'][i-1]) for i in (2,3,6,7,8)};eq=rationaldata(geo['field_equation']);fc={i:K(D['base_values'][i-1]) for i in functions}
def evalfn(coeff,u,v):
 ans=zero
 for num,den in reversed(coeff):
  d=den(v);assert d;ans=ans*u+num(v)/d
 return ans
seen=set();excluded=remaining=no_source=0;start=time.monotonic()
assert len(records)==len(F['rows'])==13513
for i,(row,saved) in enumerate(zip(records,F['rows']),1):
 assert saved[0]==i
 if row[4]==0:assert saved==[i,-1,[]];no_source+=1;continue
 b,c=elt(row[1]),elt(row[2]);t=2*b/(one-c)-one;q=((one-c)*t*t-2)/(2*t);ss=t*(q*q-one)/(t*(q*q+one)+2*q);dd=(ss+2)/q;u,v=(ss+dd)/2,(ss-dd)/2
 assert evalfn(eq,u,v)==zero
 signs={j:(1 if (evalfn(functions[j],u,v)/fc[j]).is_square() else -1) for j in functions}
 assert all(signs[j]==saved[2][j-1] for j in functions)
 passed=all(__import__('math').prod(signs.get(j+1,0)**a for j,a in enumerate(combo) if a)==1 for combo in D['annihilating_combinations'])
 assert saved[1]==int(passed)
 if passed:remaining+=1
 else:excluded+=1
 if i%2000==0:print('SAGE_CHARACTER_PROGRESS',i,round(time.monotonic()-start,2),flush=True)
assert (no_source,excluded,remaining)==(1644,11173,696)
print('NONIC48_SAGE_CHARACTER_REPLAY_PASS '+json.dumps({'records':13513,'all_function_evaluations_independently_recomputed':True,'source_free_orbits':no_source,'quadratic_exclusions':excluded,'full_queries':remaining,'complete':True}))
