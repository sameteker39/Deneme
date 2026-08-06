import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from efa_core import smooth, efa, communalities
from final_engine import msa_vec
pd.set_option('display.width',250)
d=pickle.load(open('D.pkl','rb')); R=smooth(d['R'],0.02); cols=d['cols']; p=21
m=msa_vec(d['R'])
for k in [2,3,4,5]:
    L,Phi=efa(R,k,'oblimin','minres'); h2=communalities(L,Phi)
    T=pd.DataFrame(L.round(2),index=cols,columns=['F%d'%(i+1) for i in range(k)])
    T2=T.copy()
    for c in T.columns: T2[c]=[('%.2f'%v if abs(v)>=.32 else '  .  ') for v in T[c]]
    T2['h2']=h2.round(2); T2['MSA']=m.round(2)
    print("\n===== k=%d ====="%k); print(T2.to_string())
    print("Phi:",np.round(Phi[np.triu_indices(k,1)],2))
    A=np.abs(L); pri=A.max(1)
    print("yuk<.40: %d | h2<.30: %d | faktor buyuklukleri: %s"%(
        (pri<.40).sum(),(h2<.30).sum(), pd.Series(A.argmax(1)).value_counts().sort_index().to_dict()))
