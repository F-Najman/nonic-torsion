#!/usr/bin/env python3
"""Independent read-only Sage/M4RI ranks of retained binary matrices."""
from pathlib import Path
import gzip,json,time
from sage.all import GF,matrix
H=Path(__file__).resolve().parents[1]/'data'/'quotient26';t=time.monotonic();cases=json.loads(gzip.open(H/'canonical_ring_matrices.json.gz','rt').read());assert len(cases)==60
for c in cases:
 r=[int(a,16)for a in c['rows']];n=c['columns'];assert len(r)==75 and n in(128,132);assert matrix(GF(2),[[a>>j&1 for j in range(n)]for a in r]).rank()==75
print('NONIC71_M4RI_READONLY_PASS cases60 rank75 elapsed',time.monotonic()-t)
