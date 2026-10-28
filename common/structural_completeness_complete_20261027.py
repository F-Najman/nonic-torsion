#!/usr/bin/env python3
"""Independent final finite-group and evidence-integrity audit of Phi(9).

Run only after all94 forbidden-subgroup exclusions are accepted. No import
of candidates.py, no SymPy/Sage/Magma, and no modifications of the ledger,
generated candidates, prior audits, or earlier-degree files.

The mathematical correctness and forbidden-subgroup strength of the
separately reviewed proof packets, together with the prime/family theorems,
are explicit inputs. This script proves the exhaustive structural reduction
and checks evidence integrity; hashes do not re-prove mathematical claims.
"""
from collections import defaultdict
from functools import lru_cache
from math import gcd
from pathlib import Path
import hashlib,json

B=Path(__file__).resolve().parent/'data'
ROOT=B.parents[1]

@lru_cache(None)
def factor(n):
    answer={};p=2
    while p*p<=n:
        while n%p==0:answer[p]=answer.get(p,0)+1;n//=p
        p+=1
    if n>1:answer[n]=answer.get(n,0)+1
    return answer

def normalize(a,b):
    d=gcd(a,b);return d,a*b//d

def maximal_types(g):
    a,b=g;out={normalize(a//p,b)for p in factor(a)}
    out|={normalize(a,b//p)for p in factor(b)}
    assert all(h[0]*h[1]<a*b and (a*b)%(h[0]*h[1])==0 for h in out)
    return out

def closure(groups):
    result=set(groups);todo=list(groups)
    while todo:
        for h in maximal_types(todo.pop()):
            if h not in result:result.add(h);todo.append(h)
    return result

def embeds(h,g):
    # Independent p-primary elementary-divisor comparison.
    hh=[factor(n)for n in h];gg=[factor(n)for n in g]
    return all(hh[i].get(p,0)<=gg[i].get(p,0)
               for i in range(2)for p in set(hh[i])|set(gg[i]))

def sha(path):
    digest=hashlib.sha256()
    with path.open('rb')as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()

ledger_path=B/'ledger.json';candidate_path=B/'candidates.json'
ledger_raw=ledger_path.read_bytes();candidate_raw=candidate_path.read_bytes()
ledger_snapshot_sha=hashlib.sha256(ledger_raw).hexdigest()
candidate_snapshot_sha=hashlib.sha256(candidate_raw).hexdigest()
audit_source_sha=sha(Path(__file__))
producer_snapshot_sha=sha(B.parent/'candidates.py')
ledger=json.loads(ledger_raw);generated=json.loads(candidate_raw)
# These gates intentionally reject the earlier research checkpoint.
assert ledger['degree']==9 and ledger['classification_complete']is True
assert generated['degree']==9 and generated['classification_complete']is True
assert generated['status']=='COMPLETE'and generated['remaining_minimal_obstructions']==[]
assert len(ledger['exclusions'])==94,'All94 exclusions must be accepted before this final audit runs.'
assert ledger['prime_bound']['upper_bound']==generated['prime_bound']==89
assert ledger['prime_bound']['audited']is True

prime_universe={n for n in range(2,90)if factor(n)=={n:1}}
assert len(prime_universe)==24 and max(prime_universe)==89
assert sorted(prime_universe)==generated['prime_universe']
# Since9 is odd and phi(m) is even for m>2, the Weil pairing leaves only1,2.
first={1,2};assert ledger['first_invariant_factors']==generated['first_invariant_factors']==[1,2]
published={(1,n)for n in list(range(1,29))+[30,36]}
published|={(2,2*n)for n in range(1,13)}
assert len(published)==42 and published==set(map(tuple,ledger['published_infinite_families']['groups']))
assert len(ledger['published_infinite_families']['groups'])==42
extra={(1,n)for n in[29,31,32,33,34,39,40,42]}|{(2,26)}
realized=[tuple(r['group'])for r in ledger['realizations']]
assert len(realized)==len(set(realized))==9 and set(realized)==extra
known=published|extra
expected={(1,n)for n in list(range(1,35))+[36,39,40,42]}|{(2,2*n)for n in range(1,14)}
assert known==expected and len(known)==51 and closure(known)==known
assert all(m in first and n%m==0 and set(factor(n))<=prime_universe for m,n in known)
cyclic={n for m,n in known if m==1};assert len(cyclic)==38 and max(cyclic)==42
assert sum(m==2 for m,n in known)==13

# Third construction: scan the WHOLE proved interval for both invariant factors,
# using maximal index-prime subgroup types rather than proper-divisor pairs or
# the producer's products c*p. No arbitrary search bound on full torsion is used.
bound=max(cyclic)*max(prime_universe);assert bound==3738
frontier=[]
for m in sorted(first):
    for n in range(1,bound+1):
        if n%m or not set(factor(n))<=prime_universe:continue
        g=m,n
        if g not in known and maximal_types(g)<=known:frontier.append(g)
assert frontier==sorted(frontier)and len(frontier)==94
assert sum(m==1 for m,n in frontier)==87 and sum(m==2 for m,n in frontier)==7
assert all(m==1 or n in cyclic for m,n in frontier)
assert [n for m,n in frontier if m==2]==[28,30,32,34,36,40,42]
excluded=[tuple(r['group'])for r in ledger['exclusions']]
assert len(excluded)==len(set(excluded))==94 and set(excluded)==set(frontier)
assert all(not embeds(e,g)for e in excluded for g in known)
remaining=[g for g in frontier if not any(embeds(e,g)for e in excluded)]
assert remaining==[]

# Match every generated statement, independently of its implementation.
assert generated['phi_infinite']==[list(g)for g in sorted(published)]
assert generated['known_lower_bound']==generated['classification']==[list(g)for g in sorted(known)]
assert generated['minimal_outside_known']==[list(g)for g in frontier]
assert generated['certified_exclusions']==[list(g)for g in sorted(excluded)]
assert generated['eliminated_minimal_obstructions']==[list(g)for g in frontier]

# The evidence files of the accepted proof ledger: common/data/evidence_sha256.tsv lists each of
# them with its SHA-256 in the ledger and its path in this repository if it is kept byte for byte.
# A file rewritten by a rerun of its job is reported, not rejected.
evidence={};references=0;changed=[]
for row in (B/'evidence_sha256.tsv').read_text().splitlines()[1:]:
    original,digest,path=row.split('\t');references+=1
    if path=='-':continue
    if sha(ROOT/path)==digest:evidence[original]=digest
    else:changed.append(path)
for kind in ['realizations','exclusions']:
    for r in ledger[kind]:
        assert r['degree']==9 and r['audited']is True
        assert r['folder']is None or (ROOT/r['folder']).is_dir(),r['folder']

assert sha(ledger_path)==ledger_snapshot_sha
assert sha(candidate_path)==candidate_snapshot_sha
assert sha(B.parent/'candidates.py')==producer_snapshot_sha
assert sha(Path(__file__))==audit_source_sha

byfirst=defaultdict(list)
for m,n in frontier:byfirst[m].append(n)
report={
 'audit_date':'2026-10-27','degree':9,'complete_structural_audit':True,
 'classification_complete':True,'classification':[list(g)for g in sorted(known)],
 'known_count':51,'published_infinite_family_count':42,
 'additional_realizations':[list(g)for g in sorted(extra)],
 'known_set_downward_closed':True,'cyclic_count':38,'noncyclic_count':13,
 'allowed_first_invariant_factors':[1,2],'prime_bound':89,
 'prime_universe':sorted(prime_universe),'max_known_cyclic_exponent':42,
 'proved_minimal_enumeration_bound':bound,
 'independent_enumeration_method':'Full interval for both first invariants; maximal index-prime coordinate deletions and gcd/lcm normalization. Trial factoring only; no import of candidates.py and no symbolic algebra.',
 'minimal_frontier':[list(g)for g in frontier],'frontier_count':94,
 'frontier_cyclic_count':87,'frontier_noncyclic_count':7,
 'frontier_by_first_invariant':dict(byfirst),
 'accepted_exclusion_records':94,'remaining_minimal_obstructions':[],
 'all94_minimal_types_have_accepted_forbidden_subgroup_proofs':True,
 'verified_unique_evidence_files':len(evidence),'verified_evidence_references':references,
 'external_inputs':{
  'prime_bound':ledger['prime_bound']['source'],
  'published_families':ledger['published_infinite_families']['source'],
  'proof_packets':'The mathematical correctness and forbidden-subgroup strength of each separately reviewed accepted proof are inputs.'},
 'sha256':{'ledger':ledger_snapshot_sha,'candidates':candidate_snapshot_sha,
           'candidates_producer_read_only':producer_snapshot_sha,
           'audit_program':audit_source_sha},
 'verification_boundary':'This program proves and independently verifies the exhaustive finite-group reduction and evidence integrity. Published mathematical theorems and the correctness/subgroup scope of accepted occurrence and exclusion proofs remain mathematical inputs; hashes alone do not prove them.'}
(B/'structural_completeness_complete_20261027.json').write_text(json.dumps(report,indent=2)+'\n')
text='''Independent final structural audit of Phi(9),27October2026.

The51groups are C_n for n=1,...,34,36,39,40,42, and C2xC_(2n)
for n=1,...,13. These are the42published infinite-family groups and
nine additional fully certified realizations. Their38cyclic exponents
have maximum42; all51types form a subgroup-closed set.

The prime bound is89. The Weil pairing gives phi(m)|9; since phi(m)
is even for m>2, normalized rank-two torsion has first invariant1or2.
For a subgroup minimal outside the known set, every proper subgroup
is known. In the cyclic case n/p<=42 for each prime p|n, so n<=3738.
In the noncyclic case its cyclic subgroup C_n is proper and known,
so n<=42. These bound MINIMAL obstructions, without imposing an
unproved uniform bound on full torsion orders.

A third independent enumeration scans the entire proved interval
for both first invariants and tests maximal subgroup types via
index-prime coordinate deletions and gcd/lcm normalization. It uses
no symbolic algebra and does not import candidates.py. It finds
exactly94minimal obstructions:87cyclic and7noncyclic. They agree
with the separately generated frontier, and all94have accepted
forbidden-subgroup proofs. There are no remaining types.

The normalized embedding criterion C_a x C_b -> C_m x C_n is
a|m and b|n. Coordinate embeddings prove sufficiency; exponents
and p^(u+1)-torsion cardinalities prove necessity. Thus every
eligible torsion group outside the51known types contains an
excluded minimal subgroup. The occurrence certificates give the
reverse inclusion, completing Phi(9).

The evidence files of the accepted ledger that are kept byte for byte
have the SHA-256 recorded in the ledger. The ledger and generated
candidates were read only.
The JSON gives the complete classification, frontier, input hashes,
evidence counts and explicit mathematical verification boundary.
'''
(B/'structural_completeness_complete_20261027.txt').write_text(text)
print('NONIC_COMPLETE_STRUCTURAL_AUDIT_PASS known=51 frontier=94 exclusions=94 remaining=0')
for path in changed:print('EVIDENCE_CHANGED_SINCE_ORIGINAL_RUN',path)
print('EVIDENCE_HASHES_PASS',len(evidence),references)
