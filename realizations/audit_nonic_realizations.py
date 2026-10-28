#!/usr/bin/env python3
"""Strict independent input/coverage/gcd/hash audit for nonic realizations."""
import hashlib,json,re
from fractions import Fraction
from math import gcd
from pathlib import Path
import sympy as s
H=Path(__file__).resolve().parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
data=json.loads((H/'common/data/inputs/vanhoeij_nonic_input.json').read_text())
assert sha(H/'common/data/inputs/LowDegreePlaces')==data['source_sha256']
examples=data['examples'];expected=[29,31]+[32]*8+[33]+[34]*4+[39,40,42,42]
assert [e['N'] for e in examples]==expected
source=(H/'common/data/inputs/LowDegreePlaces').read_text().splitlines()
t,x,y=s.symbols('t x y')
def poly(cs):return sum(s.Rational(c)*t**i for i,c in enumerate(cs))
for i,e in enumerate(examples,1):
 assert e['index']==i and source[e['source_line']-1]==e['source_text']
 f=poly(e['f']);xp=poly(e['x']);yp=poly(e['y'])
 assert s.degree(f,t)==9 and s.Poly(f,t).is_irreducible
 for eq in re.search(r'\[(.*) = 0, (.*) = 0\]',e['source_text']).groups():
  z=s.sympify(eq.replace('^','**')).subs({x:xp,y:yp})
  assert s.rem(z,f,t)==0
log=(H/'realizations/logs/verify_nonic_examples.log').read_text();stdout=(H/'realizations/logs/verify_nonic_examples_stdout.log').read_text()
for text in [log,stdout]:
 assert not re.search('Runtime error|User error|Assertion failed',text)
 assert 'NONIC_VANHOEIJ_EXAMPLES_DONE count=19 complete=true' in text
full=re.findall(r'NONIC_FULL_TORSION_CERTIFIED index=(\d+) invariants=\[ (\d+) \] exact=true',log)
assert full==[(str(i),str(n)) for i,n in enumerate(expected,1)]
records=re.findall(r'NONIC_GOOD_REDUCTION index=(\d+) p=(\d+) residue_degree=(\d+) order=(\d+) gcd=(\d+)',log)
bounds={i:0 for i in range(1,20)}
for z in records:
 i,p,d,n,g=map(int,z);assert s.isprime(p) and p>2 and 1<=d<=9 and n%expected[i-1]==0
 bounds[i]=gcd(bounds[i],n);assert bounds[i]==g
 disc=s.discriminant(poly(examples[i-1]['f']),t)
 assert int(s.numer(disc))%p and int(s.denom(disc))%p
 assert all(Fraction(c).denominator%p for k in ('f','x','y') for c in examples[i-1][k])
assert [bounds[i] for i in range(1,20)]==expected
selected=json.loads((H/'realizations/data/selected_reductions.json').read_text())
direct=(H/'realizations/logs/replay_direct_counts.log').read_text()
assert not re.search('Runtime error|User error|Assertion failed',direct)
assert 'NONIC_DIRECT_COUNTS_COMPLETED examples=19 reductions=26' in direct
assert len(re.findall('NONIC_DIRECT_COUNT index=',direct))==26
for e in selected:
 b=0
 for p,d,n in e['reductions']:
  assert any(tuple(map(int,z[:4]))==(e['index'],p,d,n) for z in records)
  match=re.search(rf'NONIC_DIRECT_COUNT index={e["index"]} p={p} residue_degree={d} expected={n} observed=\[([^]]+)\] exact=true',direct)
  assert match and n in [int(x) for x in match[1].replace(',',' ').split()]
  b=gcd(b,n)
 assert b==e['N']==expected[e['index']-1]
isog=(H/'realizations/logs/verify_nonic_two_isogenies.log').read_text()
assert not re.search('Runtime error|User error|Assertion failed',isog)
assert 'NONIC_TWO_ISOGENIES_DONE targets=15 complete=true' in isog
assert isog.count('NONIC_TWO_ISOGENY_CERTIFIED')==15
# (the Codex packet manifest realizations_manifest.json, with the SHA-256 of the packet files, is not written here)
print('NONIC_REALIZATIONS_AUDIT_PASS examples=19 new_groups=8 independent_direct_reductions=26 complete=true')
