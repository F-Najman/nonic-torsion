#!/usr/bin/env python3
"""Read-only independent M4RI rank replay of explicit small certificates."""
from pathlib import Path
import gzip,json,time
from sage.all import GF,matrix
H=Path(__file__).resolve().parents[2]/'data/P89/quotient27';t0=time.monotonic()
with gzip.open(H/'canonical_ring_small_matrices.json.gz','rt')as fp:cases=json.load(fp)
assert len(cases)==96
F=GF(2)
for c in cases:
 rows=[int(x,16)for x in c['rows']];assert len(rows)==130 and c['columns']==264
 A=matrix(F,[[x>>j&1 for j in range(264)]for x in rows]);assert A.rank()==130
print('NONIC89_SMALL_MATRICES_SAGE_M4RI_PASS cases=96 rank=130 time=',time.monotonic()-t0)
