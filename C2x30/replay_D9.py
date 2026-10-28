"""Independent integer and finite-field replay of the nonic C2x30 packet.

No Magma library is called. The proved native RR and complete factorization
stages are reproducible separately; this checks their exported arithmetic,
all series classification, local lift combinations, and actual level maps.
"""
from pathlib import Path
from collections import Counter,deque
import runpy,json,itertools,sys
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
from finite_arithmetic import read_magma_data,Field,Curve,evaluate

P=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent/'logs'/'run'

def key_mod(v,iv):
    assert len(v)==len(iv)+1 and v[-1]==0
    return tuple(x%m for x,m in zip(v,iv))

def series():
    a=runpy.run_path(str(P/'classes256.py'));assert a['complete']
    R=a['relations'];n=len(R);iv=[R[i][i] for i in range(n)]
    assert all(row==[iv[i] if j==i else 0 for j in range(n+1)] for i,row in enumerate(R))
    hs=a['cusp_block_classes'];ds=a['cusp_block_degrees'];assert ds==[1,1,4,2,4,1,2,1]
    H={(0,)*n};q=deque(H)
    while q:
        x=q.popleft()
        for h in hs:
            y=tuple((x[i]+h[i])%iv[i] for i in range(n))
            if y not in H:H.add(y);q.append(y)
    assert len(H)==256
    best={};count=0
    def enum(j,left,v):
        nonlocal count
        if j==8:
            if left:return
            count+=1;k=tuple(sum(v[t]*hs[t][i] for t in range(8))%iv[i] for i in range(n))
            b=sum((v[t]//2)*ds[t] for t in range(8));h=max(1,b)
            if k not in best or best[k][0]<h:best[k]=(h,list(v),b)
            return
        for m in range(left//ds[j]+1):enum(j+1,left-m*ds[j],v+[m])
    enum(0,9,[]);assert count==956 and len(best)==256
    native=a['cases'];assert len(native)==256 and {key_mod(r[1],iv) for r in native}==H
    routes=Counter();exceptions=[]
    for idx,g,known,v,h in native:
        k=key_mod(g,iv);bh,bv,b=best[k]
        assert known and all(x>=0 for x in v) and sum(x*d for x,d in zip(v,ds))==9
        assert tuple(sum(v[t]*hs[t][i] for t in range(8))%iv[i] for i in range(n))==k
        assert bh<=h
        if h==1:route='unique_cuspidal'
        elif h==bh:
            route='pullback_with_fixed_cusps';assert b>0 and sum((x%2)*d for x,d in zip(bv,ds))==9-2*b>0
        else:
            route='exceptional_series';assert h==2 and bh==1;exceptions.append([idx,bv,h])
        routes[route]+=1
    assert routes==dict(unique_cuspidal=84,pullback_with_fixed_cusps=144,exceptional_series=28)
    inputs=read_magma_data(P/'exceptional_pencil_inputs.m')
    assert inputs['C230ExpectedCuspSignatures']==a['cusp_block_signatures']
    assert inputs['C230ExpectedCuspDegrees']==ds and inputs['C230ExceptionalCases']==exceptions
    print('C230_D9_SERIES_PASS classes=256 effective_divisors=956 unique=84 fixed_cusp=144 pencils=28')
    return exceptions

class Quotient:
    def __init__(self,R,hs):
        self.n=len(R[0])-1
        assert len(R)==self.n and all(r[-1]==0 for r in R+hs)
        self.H=hermite_normal_form(Matrix([r[:-1] for r in R+hs]).T)
        assert self.H.shape==(self.n,self.n)
        self.index=abs(int(self.H.det()))
    def key(self,v):
        assert len(v)==self.n+1 and v[-1]==0
        w=list(v[:-1])
        for j in reversed(range(self.n)):
            d=int(self.H[j,j]);assert d>0;q=w[j]//d
            for i in range(j+1):w[i]-=q*int(self.H[i,j])
        return tuple(w)


def fibres_and_maps(exceptions,target_folder=None):
    a=read_magma_data(P/'exceptional_fibres_data.m');assert a['C230FibresComplete']
    ps=a['C230Points'];cs=a['C230Candidates'];fs=a['C230Fibres']
    assert len(ps)==78 and len(cs)==42 and len(fs)==224
    assert len({repr(p) for p in ps})==78 and len({repr(c) for c in cs})==42
    assert {(r[0],r[1]) for r in fs}=={(e[0],i) for e in exceptions for i in range(1,9)}
    for idx,sid,nc,valid,parts,exported in fs:
        if not valid:
            assert not exported and any(dd==[3] and m%3 for de,m,dd,ids,n in parts)
            continue
        if nc==0:assert not exported;continue
        current=[[]]
        for de,m,dd,ids,n in parts:
            if dd==[1,1,1]:
                assert len(ids)==3 and len(set(ids))==3 and all(ps[i-1][0]==de for i in ids)
                opts=[[[ids[i],v] for i,v in enumerate([aa,bb,m-aa-bb]) if v] for aa in range(m+1) for bb in range(m-aa+1)]
            else:
                assert dd==[3] and m%3==0 and len(ids)==3 and len(set(ids))==1 and ps[ids[0]-1][0]==3*de
                opts=[[[ids[0],m//3]]]
            assert len(opts)==n
            current=[old+opt for old in current for opt in opts]
        actual=set()
        for cand in current:
            combined=Counter()
            for i,m in cand:combined[i]+=m
            clean=[[i,combined[i]] for i in sorted(combined)]
            assert sum(ps[i-1][0]*m for i,m in clean)==nc
            actual.add(cs.index(clean)+1)
        assert set(exported)==actual
    assert {i for row in fs for i in row[-1]}==set(range(1,43))
    assert Counter(sum(ps[i-1][0]*m for i,m in c) for c in cs)=={8:18,9:24}
    assert all(len(c)==1 and c[0][1]==1 for c in cs)
    active=sorted({r[0] for c in cs for r in c});assert len(active)==42
    b=runpy.run_path(str((target_folder or P)/'three_quotients.py'));assert b['complete'] and not b['remaining']
    targets={r[0]:r[1] for r in b['targets']};assert len(targets)==21
    assert len(b['source_maps'])==42 and len(b['candidate_tests'])==42
    fieldcurves={}
    for pid,(de,bits,pol,coeff) in enumerate(ps,1):
        K=Field(pol);assert K.n==de and de in [2,4,8,9]
        bc=[evaluate(v,K.z) for v in coeff];assert bits[0]*bc[0]+bits[1]*bc[1]==K.z
        assert evaluate(pol,K.z)==0
        E=Curve(K,*bc);P2=(K.zero,K.zero);Q=(bc[0],bc[0]);E.check(P2);E.check(Q)
        assert E.mul(2,P2) is None and E.mul(30,Q) is None
        assert all(E.mul(30//p,Q) is not None for p in [2,3,5]) and P2!=E.mul(15,Q)
        fieldcurves[pid]=(K,E,P2,Q)
    tests=0
    for pos,pid in enumerate(active):
        K,E,P2,Q=fieldcurves[pid];de=K.n
        variants=[Q,E.add(P2,E.mul(16,Q)),E.add(P2,Q)]
        maps=b['source_maps'][pos];assert [r[0] for r in maps]==[1,2,3]
        for (variant,ti,passed),R in zip(maps,variants):
            assert E.mul(30,R) is None and all(E.mul(30//p,R) is not None for p in [2,3,5])
            pair=E.tate_xy(R);desc=targets[ti];e=desc[0];assert de%e==0 and len(desc)==4+3*e
            bits=desc[1:3];pol=desc[3:4+e];coeff=[desc[4+e:4+2*e],desc[4+2*e:4+3*e]]
            gen=bits[0]*pair[0]+bits[1]*pair[1];assert evaluate(pol,gen)==0
            assert all(gen**(7**d)!=gen for d in range(1,e) if e%d==0)
            assert [evaluate(v,gen) for v in coeff]==pair
            tests+=1
    assert tests==126
    print('C230_D9_FIBRES_PASS members=224 points=78 source_candidates=42 nc8=18 nc9=24')
    print('C230_D9_FINITE_ARITHMETIC_PASS primitive_levels=78 quotient_maps=126 targets=21')
    Q=Quotient(b['relations'],b['rational_cusp_generators']);assert Q.index==4352
    zero=(0,)*Q.n
    allowed={0:{zero},1:{Q.key(v) for d,v in zip(b['finite_cusp_degrees'],b['finite_cusp_classes']) if d==1}}
    assert sorted(b['finite_cusp_degrees'])==[1]*16+[4]*4
    targetdata={}
    for ti,desc,v,passed in b['targets']:
        de=desc[0];assert de in [8,9]
        assert bool(passed)==(Q.key([-x for x in v]) in allowed[9-de])
        targetdata[ti]=bool(passed)
    pointpassed={}
    for pid,maps in zip(active,b['source_maps']):
        assert all(bool(passed)==targetdata[ti] for _,ti,passed in maps)
        pointpassed[pid]=all(targetdata[ti] for _,ti,_ in maps)
    seen=set()
    for idx,pid,passed in b['candidate_tests']:
        assert 1<=idx<=42 and idx not in seen;seen.add(idx)
        assert pid==cs[idx-1][0][0] and bool(passed)==pointpassed[pid] and not passed
    assert seen==set(range(1,43))
    print('C230_D9_INTEGER_PASS quotient_index=4352 finite_cusp_shifts1='+str(len(allowed[1]))+' candidates=42 remaining=0')

if __name__=='__main__':
    exc=series();fibres_and_maps(exc)
    control=P.parent/'control'
    if P.name=='run' and (control/'three_quotients.py').is_file():
        fibres_and_maps(exc,control)
        print('C230_D9_ORDINARY_NATIVE_FORWARD_MAP_CONTROL_PASS')
    print('C230_D9_INDEPENDENT_REPLAY_PASS complete=true')
