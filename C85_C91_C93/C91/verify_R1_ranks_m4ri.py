#!/usr/bin/env python3
"""Independent Sage/M4RI ranks, without native launch or output mutation."""
from pathlib import Path
import gzip,json,time
from sage.all import GF,matrix
H=Path(__file__).resolve().parents[1]/'data/C91';t=time.monotonic();rr=json.loads(gzip.open(H/'canonical_R1_matrices.json.gz','rt').read());assert len(rr)==288
for c in rr:
 rows=[int(v,16)for v in c['rows']];n=c['columns'];assert len(rows)==265 and n>=527 and all(v.bit_length()<=n for v in rows);assert matrix(GF(2),[[v>>j&1 for j in range(n)]for v in rows]).rank()==265
print('NONIC91_R1_M4RI_RANK_PASS cases288 rank265 elapsed',time.monotonic()-t)
