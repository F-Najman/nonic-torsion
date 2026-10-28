"""Independent full66-orbit arithmetic and exact diagonal stabilizer replay.
This is a partial necessary-condition certificate, not an exclusion.
"""
from pathlib import Path
from collections import Counter
import ast,runpy
from sympy import Matrix
from sympy.matrices.normalforms import hermite_normal_form
import sys;sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'common'))
from finite_arithmetic_char3 import read_magma_data,Field,Curve,evaluate
P=Path(__file__).resolve().parent/'data'

def read_data(q):
    s=q.read_text().replace('[*','[').replace('*]',']').replace('true','True').replace('false','False')
    return ast.literal_eval(s)

def frobkey(vs):
    keys=[]
    for i in range(9):
        keys.append(tuple(z for v in vs for z in v.v));vs=[v**3 for v in vs]
    return min(keys)

def inside(v,H):
    y=H.inv()*Matrix(v[:-1]);return all(a.q==1 for a in y)

def main():
    inv,rels,hg,rs,remain,fp,inputs,labels,orbits=read_data(P/'N26_D9_threequotients_arithmetic.txt')
    assert inv==[665,1995,0] and len(rs)==66 and len(inputs)==2328 and len(set(map(tuple,inputs)))==2328
    assert len(labels)==2328 and set(labels)==set(range(1,67)) and len(orbits)==66
    assert Counter(len(set(map(tuple,o))) for o in orbits)=={36:64,12:2}
    assert [r[0] for r in remain]==[21,58]
    H=hermite_normal_form(Matrix([v[:-1] for v in rels+hg]).T);assert abs(int(H.det()))==5
    assert all(inside(v,H) for v in hg)
    native_stabs={r[0]:r for r in read_data(P/'N26_D9_stabilizers_arithmetic.txt')}
    K=Field(rs[0][1]);assert K.n==9
    assert all(r[1]==rs[0][1] for r in rs)
    targets={};source_union=set();native_inputs={tuple(x) for x in inputs}
    for rec in rs:
        idx,pol,vals,maps,passed=rec;assert idx in range(1,67)
        b,c=[evaluate(v,K.z) for v in vals];E=Curve(K,b,c);T=(K.zero,K.zero);Q=(b,b)
        E.check(T);E.check(Q);assert E.mul(2,T) is None and E.mul(26,Q) is None
        assert E.mul(13,Q) is not None and E.mul(2,Q) is not None and T!=E.mul(13,Q)
        twos=[T,E.mul(13,Q),E.add(T,E.mul(13,Q))];V=E.mul(14,Q)
        multiples=[E.mul(a,V) for a in range(1,7)]
        actual=[];stabs=[];original=frobkey([b,c])
        for pi in range(1,4):
            for ri in range(1,4):
                if pi==ri:continue
                U=twos[pi-1]
                for aa in range(1,7):
                    W=E.add(twos[ri-1],multiples[aa-1]);xx=W[0]-U[0];yy=W[1];iy=(yy*yy).inv()
                    pair=[xx**3*iy,(c+3*U[0])*xx*xx*iy];key=frobkey(pair);actual.append(key)
                    if key==original:
                        phases=[i for i in range(9) if pair==[b**(3**i),c**(3**i)]];assert len(phases)==1
                        stabs.append([pi,ri,aa,phases[0]])
        assert set(actual)==set(map(tuple,orbits[idx-1]))
        assert original in actual and len(set(actual)) in [12,36]
        assert set(actual)=={tuple(inputs[j]) for j in range(2328) if labels[j]==idx}
        source_union|=set(actual)
        if idx in [21,58]:
            native=native_stabs[idx]
            assert native[1]==pol and native[2]==vals
            assert stabs==[[v[0],v[1],v[2],v[7]] for v in native[3]]
            assert len(stabs)==3 and stabs[1][:3]==[2,3,3]
            assert stabs[1][3]==(3 if idx==21 else 6)
            U=E.mul(13,Q);W=E.add(T,E.mul(3,Q))
            U2=E.mul(13,W);W2=E.add(U,E.mul(3,W))
            U3=E.mul(13,W2);W3=E.add(U2,E.mul(3,W2))
            assert U3==T and W3==Q
        else:assert stabs==[[1,2,1,0]]
        variants=[Q,E.add(T,Q),E.add(T,E.mul(14,Q))]
        assert len(maps)==3
        flags=[]
        for variant,mp in enumerate(maps,1):
            vn,ti,e,mult,cl,flag=mp;assert vn==variant and e==9 and mult==1
            pair=E.tate_xy(variants[variant-1]);key=frobkey(pair)
            if ti in targets:assert targets[ti]==(key,cl)
            else:targets[ti]=(key,cl)
            flag2=inside(cl,H);assert bool(flag)==flag2;flags.append(flag2)
        assert bool(passed)==all(flags)==(idx in [21,58])
        print('N26_D9_PARTIAL_REPLAY_ORBIT',idx,'passes='+str(bool(passed)),flush=True)
    assert source_union==native_inputs
    # Independent source Newton identities establish exactly2328NC places ofdegree9.
    source=read_data(P/'degree8/N26_D8_class_arithmetic.txt');sf=source[5];traces=[];closed=[]
    for n in range(1,10):
        tr=-n*sf[50-n]-sum(sf[50-j]*traces[n-j-1] for j in range(1,n));traces.append(tr)
        num=3**n+1-tr-sum(d*closed[d-1] for d in range(1,n) if n%d==0);assert num%n==0;closed.append(num//n)
    assert closed==[18,0,6,0,36,122,288,936,2328]
    assert len(fp)==21 and sum(fp)==665*1995
    last=read_data(P/'N26_D9_last2_Hecke_arithmetic.txt');assert last[1]==[21,58]
    assert last[0]==[[i,q,True] for i in [21,58] for q in [5,7,11]]
    print('N26_D9_PARTIAL_REPLAY_PASS complete_finite_orbits=66 places=2328 maps=198 surviving_orbits=21,58')
    print('N26_D9_DIAGONAL_STABILIZER_PASS h(P,Q)=(13Q,P+3Q) order=3 phases=3,6 quotient_closed_degree=3')
    print('N26_D9_EXCLUSION_PROVED=false')

if __name__=='__main__':main()
