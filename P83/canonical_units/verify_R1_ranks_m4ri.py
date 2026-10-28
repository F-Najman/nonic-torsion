#!/usr/bin/env python3
"""Independent Sage/M4RI ranks, without native launch or output mutation."""
from pathlib import Path
import gzip,json,time
from sage.all import GF,matrix
H=Path(__file__).resolve().parents[1]/'data/canonical_units';t=time.monotonic();rr=json.loads(gzip.open(H/'canonical_R1_matrices.json.gz','rt').read());assert len(rr)==410
for c in rr:
 rows=[int(v,16)for v in c['rows']];n=c['columns'];assert len(rows)==247 and n>=491 and all(v.bit_length()<=n for v in rows);assert matrix(GF(2),[[v>>j&1 for j in range(n)]for v in rows]).rank()==247
print('NONIC83_R1_M4RI_RANK_PASS cases410 rank247 elapsed',time.monotonic()-t)
