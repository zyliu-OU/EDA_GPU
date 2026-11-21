import json, glob
from pathlib import Path
import numpy as np
def rx(p):
    L=Path(p).read_text().split('\n'); n=int(L[0]); s=[];c=[]
    for l in L[2:2+n]: t=l.split(); s.append(t[0]); c.append([float(x) for x in t[1:4]])
    return s,np.array(c)
cases=json.load(open('reports/hplc_patch_extended_2026-11-11/cases.json'))
def geom(X,s):
    C1=X[[0,26,52]]; cen=C1.mean(0); nrm=np.cross(C1[1]-C1[0],C1[2]-C1[0]); nrm/=np.linalg.norm(nrm)
    if np.dot(X[[7,33,59]].mean(0)-cen,nrm)<0: nrm=-nrm
    out={}
    # per carbon level k: mean pair distance between chains, mean radial distance from axis
    pair=[];rad=[]
    for k in range(8):
        P=X[[k,26+k,52+k]]
        pair.append(np.mean([np.linalg.norm(P[i]-P[j]) for i,j in ((0,1),(0,2),(1,2))]))
        v=P-cen; rad.append(np.mean(np.linalg.norm(v-np.outer(v@nrm,nrm),axis=1)))
    out['pair']=pair; out['rad']=rad
    # per-chain radial of C8 and max pair at top half
    P=X[[7,33,59]]; out['pairs8']=sorted(np.linalg.norm(P[i]-P[j]) for i,j in ((0,1),(0,2),(1,2)))
    out['top']=np.mean(pair[4:])
    # hollow free radius: min distance of any patch atom (C5-C8 region height) to axis  -> use all patch atoms with height>3
    v=X[:78]-cen; hh=v@nrm; lat=np.linalg.norm(v-np.outer(hh,nrm),axis=1)
    out['free']=lat[hh>3].min()
    # chain-chain min H..H / closest contacts between different chains
    ch=[list(range(i*26,(i+1)*26)) for i in range(3)]
    D=np.linalg.norm(X[:78,None]-X[None,:78],axis=-1)
    out['ccmin']=[D[np.ix_(ch[i],ch[j])].min() for i,j in ((0,1),(0,2),(1,2))]
    # guest heavy atoms inside triangle prism (lat<3, below C8 mean height)
    if len(s)>78:
        g=[i for i in range(78,len(s)) if s[i]!='H']; v=X[g]-cen; gh=v@nrm; gl=np.linalg.norm(v-np.outer(gh,nrm),axis=1)
        top=((X[[7,33,59]]-cen)@nrm).mean()
        out['gin']=int(((gl<3.0)&(gh<top)).sum()); out['gheavy']=len(g)
    # gauche count in chains
    def dih(a,b,c,d):
        b0=a-b;b1=c-b;b2=d-c;b1/=np.linalg.norm(b1);v=b0-(b0@b1)*b1;w=b2-(b2@b1)*b1
        return np.degrees(np.arctan2(np.cross(b1,v)@w,v@w))
    out['gauche']=sum(abs(dih(*X[[o+k,o+k+1,o+k+2,o+k+3]]))<120 for o in (0,26,52) for k in range(5))
    return out
s0,X0=rx('examples/hplc/c8x3_patch.xyz'); g0=geom(X0,s0)
print('SEED pair per C level', np.round(g0['pair'],2), 'free', round(g0['free'],2))
rows={}
for d in sorted(glob.glob('reports/hplc_patch_*/patch_*_[0-9]*')):
    name=Path(d).name
    if name in rows or name not in cases or 'augdz' in d: continue
    s,X=rx(f'{d}/optimized.xyz'); g=geom(X,s); c=cases[name]
    g.update(r=c['r'],h=c['h'],tot=c['tot'],ins=(c['r']<2.5 and c['h']<8),dE=c.get('dE'))
    rows[name]=g
print(len(rows),'cases')
print(f"{'case':30s} ins  r    h    tot   top   C8pairs(min,mid,max)   free  gin gauche  pair C3..C8")
for k,g in sorted(rows.items(), key=lambda kv:(not kv[1]['ins'],kv[0])):
    print(f"{k[6:]:30s} {'Y' if g['ins'] else '.'} {g['r']:4.1f} {g['h']:4.1f} {g['tot']:6.2f} {g['top']:5.2f}  {g['pairs8'][0]:5.2f} {g['pairs8'][1]:5.2f} {g['pairs8'][2]:5.2f}  {g['free']:5.2f} {g['gin']:2d}/{g['gheavy']:<2d} {g['gauche']:2d}   "+' '.join(f'{x:4.1f}' for x in g['pair'][2:]))
import statistics as st
for lab,sel in (('inserted',True),('other',False)):
    R=[g for g in rows.values() if g['ins']==sel]
    print(lab,len(R),'top-half mean pair %.2f +- %.2f'%(st.mean(g['top'] for g in R),st.pstdev(g['top'] for g in R)),'C8 mean %.2f'%st.mean(np.mean(g['pairs8']) for g in R),'C8 max-pair mean %.2f'%st.mean(g['pairs8'][2] for g in R),'free %.2f'%st.mean(g['free'] for g in R),'gauche %.1f'%st.mean(g['gauche'] for g in R))
json.dump({k:{a:(list(map(float,b)) if isinstance(b,(list,np.ndarray)) else (bool(b) if isinstance(b,(bool,np.bool_)) else float(b))) for a,b in g.items()} for k,g in rows.items()}, open('reports/hplc_patch_splay_2026-11-16/splay.json','w'), indent=1)
