import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
from efa_core import smooth
d=pickle.load(open('K17.pkl','rb')); dfc=d['dfc']; K17=d['keep']
def fa_eig(R):
    Ri=np.linalg.inv(smooth(R,0.02)); smc=1-1/np.diag(Ri)
    Rr=R.copy(); np.fill_diagonal(Rr,smc); return np.linalg.eigvalsh(Rr)[::-1]
rng=np.random.default_rng(808); B=400
for nm,cs in [('17 madde',K17),('15 madde (-M16,-M21)',[c for c in K17 if c not in('M16','M21')])]:
    sub=dfc[cs]; p=len(cs); V=sub.values
    Ro,_=polychoric_matrix(sub); Ro=smooth(Ro,0.02)
    opc=np.linalg.eigvalsh(Ro)[::-1]; ofa=fa_eig(Ro)
    pc=np.zeros((B,p)); fa=np.zeros((B,p))
    for b in range(B):
        P=V.copy()
        for j in range(p):
            col=P[:,j]; ok=~np.isnan(col); col[ok]=rng.permutation(col[ok]); P[:,j]=col
        Rb,_=polychoric_matrix(pd.DataFrame(P,columns=cs)); Rb=smooth(Rb,0.02)
        pc[b]=np.linalg.eigvalsh(Rb)[::-1]; fa[b]=fa_eig(Rb)
    print("\n=== Paralel Analiz: %s (B=%d) ==="%(nm,B))
    print(pd.DataFrame({'k':range(1,6),'Gozl_PC':opc[:5].round(3),'Rast95_PC':np.percentile(pc,95,0)[:5].round(3),
                        'Gozl_FA':ofa[:5].round(3),'Rast95_FA':np.percentile(fa,95,0)[:5].round(3)}).to_string(index=False))
    for lbl,o,s in [('PC-95%',opc,np.percentile(pc,95,0)),('FA-95%',ofa,np.percentile(fa,95,0))]:
        print("   %s -> k=%d"%(lbl,int(np.argmax(o<=s)) if (o<=s).any() else p))
