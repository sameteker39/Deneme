import pickle, numpy as np, pandas as pd, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from efa_core import smooth, efa, fit_indices
d=pickle.load(open('R.pkl','rb')); R=smooth(d['R']); cols=d['cols']; N=306
rows=[]
for k in range(1,9):
    L,Phi=efa(R,k,'oblimin','minres')
    fi=fit_indices(R,L,k,N,Phi)
    h2=np.diag(L@Phi@L.T)
    prim=np.abs(L).max(1)
    rows.append(dict(k=k, **{a:fi[a] for a in ['chi2','df','p','RMSEA','CFI','TLI','SRMR']},
                     dusuk_yuk=int((prim<.40).sum()), dusuk_h2=int((h2<.30).sum()),
                     acik_var=round(100*(L**2).sum()/31,1)))
t=pd.DataFrame(rows)
t['chi2']=t['chi2'].round(1); t['p']=t['p'].round(4)
for c in ['RMSEA','CFI','TLI','SRMR']: t[c]=t[c].round(3)
print("=== 31 MADDE, oblimin dondurmeli MINRES/ULS (polikorik) ===")
print(t.to_string(index=False))
