import numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from efa_core import smooth, efa, fit_indices, communalities

CRIT = dict(h2=.30, load=.40, cross_sec=.32, cross_diff=.20, msa=.60, minitems=3)

def msa_vec(R):
    Ri=np.linalg.inv(smooth(R,0.02)); D=np.sqrt(np.diag(Ri))
    P=-Ri/np.outer(D,D); np.fill_diagonal(P,1)
    r2=(R**2).sum(0)-1; p2=(P**2).sum(0)-1
    return r2/(r2+p2)

def solution(R, cols, k, N):
    Rs=smooth(R,0.02); L,Phi=efa(Rs,k,'oblimin','minres')
    h2=communalities(L,Phi); A=np.abs(L); o=np.argsort(-A,1)
    n=len(cols); idx=np.arange(n)
    prim=A[idx,o[:,0]]; sec=A[idx,o[:,1]] if k>1 else np.zeros(n)
    d=pd.DataFrame({'Madde':cols,'Faktor':o[:,0]+1,'Yuk':prim,'2.Yuk':sec,
                    'Fark':prim-sec,'h2':h2,'MSA':msa_vec(R)})
    d['Yuk_isaretli']=L[idx,o[:,0]]
    return d, L, Phi, fit_indices(Rs,L,k,N,Phi)

def violations(d, crit=CRIT):
    v=[]
    for _,r in d.iterrows():
        if r['h2']<crit['h2']:      v.append((r['Madde'],'h2=%.2f < .30'%r['h2'],1,r['h2']))
        elif r['Yuk']<crit['load']: v.append((r['Madde'],'yuk=%.2f < .40'%r['Yuk'],2,r['Yuk']))
        elif r['2.Yuk']>=crit['cross_sec'] and r['Fark']<crit['cross_diff']:
                                    v.append((r['Madde'],'binisik %.2f/%.2f (fark %.2f)'%(r['Yuk'],r['2.Yuk'],r['Fark']),3,r['Fark']))
        elif r['MSA']<crit['msa']:  v.append((r['Madde'],'MSA=%.2f < .60'%r['MSA'],4,r['MSA']))
    v.sort(key=lambda t:(t[2],t[3]))
    return v

def purify_at_k(R0, cols0, k, N, crit=CRIT, maxdrop=None, verbose=True):
    R=R0.copy(); cols=list(cols0); log=[]
    maxdrop = maxdrop if maxdrop else int(.5*len(cols0))
    while True:
        d,L,Phi,fi=solution(R,cols,k,N)
        v=violations(d,crit)
        if not v: return dict(R=R,cols=cols,d=d,L=L,Phi=Phi,fi=fi,log=log,ok=True)
        if len(log)>=maxdrop: return dict(R=R,cols=cols,d=d,L=L,Phi=Phi,fi=fi,log=log,ok=False)
        drop,reason=v[0][0],v[0][1]
        j=cols.index(drop); keep=[i for i in range(len(cols)) if i!=j]
        R=R[np.ix_(keep,keep)]; cols=[cols[i] for i in keep]
        log.append((drop,reason,len(cols)))
        if verbose: print("     -%-4s  %-34s -> %d madde"%(drop,reason,len(cols)))

def ordinal_alpha(R):
    p=R.shape[0]; s=R.sum()
    return p/(p-1)*(1-p/s)

def omega(L, Phi, idx):
    """faktore ait maddeler icin McDonald omega (tek faktorlu yaklasim)."""
    l=L[idx]; h=communalities(L,Phi)[idx]
    sl=l.max(1) if l.ndim>1 else l
    return (sl.sum()**2)/((sl.sum()**2)+ (1-h).sum())
