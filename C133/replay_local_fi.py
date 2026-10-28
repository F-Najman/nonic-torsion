#!/usr/bin/env python3
"""Independent GF(2) replay of the local-completion C133 certificate."""
from pathlib import Path
from itertools import product
from functools import reduce
import json,sys
P=Path(__file__).resolve().parent

def matrix(entries,n,m):
 assert len(entries)==n*m and set(entries)<={0,1}
 return [sum(entries[i*m+j]<<j for j in range(m)) for i in range(n)]
def pivots(rows):
 pp={}
 for v in rows:
  while v:
   k=v.bit_length()-1
   if k in pp:v^=pp[k]
   else:pp[k]=v;break
 return pp
def contains(pp,v):
 while v:
  k=v.bit_length()-1
  if k not in pp:return False
  v^=pp[k]
 return True
def apply(v,rows):
 a=0
 while v:
  b=v&-v;a^=rows[b.bit_length()-1];v^=b
 return a
def mul(A,B):return [apply(a,B) for a in A]
def add(*args):return [reduce(int.__xor__,v,0) for v in zip(*args)]
def equalspans(A,B):
 a,b=pivots(A),pivots(B)
 return len(a)==len(b) and all(contains(a,v) for v in B)

# main() removed: it replayed fi3_local_optimized_data.json, the incomplete output of an exploratory
# run that is not part of the proof; replay_split_fi.py uses only the functions above.
