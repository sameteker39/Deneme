import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from efa_core import smooth, efa, fit_indices, communalities
pd.set_option('display.width',200)
d=pickle.load(open('R25.pkl','rb')); R=smooth(d['R'],0.02); cols=d['cols']
for k in [2,3,4]:
    L,Phi=efa(R,k,'oblimin','minres'); h2=communalities(L,Phi)
    T=pd.DataFrame(L.round(2),index=cols,columns=['F%d'%(i+1) for i in range(k)])
    T2=T.copy()
    for c in T.columns: T2[c]=[('%.2f'%v if abs(v)>=.32 else '  .  ') for v in T[c]]
    T2['h2']=h2.round(2)
    print("\n===== k=%d (25 madde, GPA oblimin) ====="%k)
    print(T2.to_string())
    print("Phi:\n",pd.DataFrame(Phi.round(2),index=T.columns,columns=T.columns).to_string())
