#!/usr/bin/env python3
"""Read-only independent M4RI rank replay of explicit small certificates."""
from pathlib import Path
import gzip,json,time
from sage.all import GF,matrix
H=Path(__file__).resolve().parents[2]/'data/P79/quotient18';t0=time.monotonic()
with gzip.open(H/'canonical_ring_matrices.json.gz','rt')as fp:cases=json.load(fp)
assert len(cases)==48
F=GF(2)
for c in cases:
 rows=[int(x,16)for x in c['rows']];assert len(rows)==119 and c['columns']==235
 A=matrix(F,[[x>>j&1 for j in range(235)]for x in rows]);assert A.rank()==119
print('NONIC79_SMALL_MATRICES_SAGE_M4RI_PASS cases=48 rank=119 time=',time.monotonic()-t0)
