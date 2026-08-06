import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
from efa_core import smooth
d=pickle.load(open('Rc.pkl','rb')); dfc=d['dfc']
def fa_eig(R):
    Ri=np.linalg.inv(smooth(R,0.02)); smc=1-1/np.diag(Ri)
    Rr=R.copy(); np.fill_diagonal(Rr,smc); return np.linalg.eigvalsh(Rr)[::-1]
rng=np.random.default_rng(77); B=300
for name,cs in [('25 madde',[c for c in dfc.columns if c not in ['M3','M4','M12','M14','M18','M20']]),
                ('22 madde',[c for c in dfc.columns if c not in ['M3','M4','M12','M14','M18','M20','M11','M17','M29']])]:
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
    print("\n=== Paralel Analiz: %s (B=%d) ==="%(name,B))
    t=pd.DataFrame({'k':range(1,7),'Gozl_PC':opc[:6].round(3),'Rast95_PC':np.percentile(pc,95,0)[:6].round(3),
                    'Gozl_FA':ofa[:6].round(3),'Rast95_FA':np.percentile(fa,95,0)[:6].round(3)})
    print(t.to_string(index=False))
    for nm,o,s in [('PC-95%',opc,np.percentile(pc,95,0)),('FA-95%',ofa,np.percentile(fa,95,0)),
                   ('FA-ort',ofa,fa.mean(0))]:
        print("   %s -> k=%d"%(nm, int(np.argmax(o<=s)) if (o<=s).any() else p))
