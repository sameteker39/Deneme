import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
from efa_core import smooth, efa, fit_indices, communalities
from robust_fit import robust_fit
from final_engine import msa_vec
pd.set_option('display.width',220)
g=pickle.load(open('gamma.pkl','rb')); boot=g['boot']; cols31=g['cols']; N=306
iu31=np.triu_indices(31,1); pair={}
for t,(i,j) in enumerate(zip(*iu31)): pair[(cols31[i],cols31[j])]=t
def gamma_for(cs):
    sel=[pair[(cs[i],cs[j])] if (cs[i],cs[j]) in pair else pair[(cs[j],cs[i])]
         for i in range(len(cs)) for j in range(i+1,len(cs))]
    return np.cov(boot[:,sel].T)*310
d=pickle.load(open('Rc.pkl','rb')); dfc=d['dfc']
C22=[c for c in cols31 if c not in ['M3','M4','M12','M14','M18','M20','M11','M17','M29']]
R22,_=polychoric_matrix(dfc[C22]); Rs=smooth(R22,0.02)
for k in [2,3]:
    L,Phi=efa(Rs,k,'oblimin','minres'); h2=communalities(L,Phi)
    T=pd.DataFrame(L.round(2),index=C22,columns=['F%d'%(i+1) for i in range(k)])
    T2=T.copy()
    for c in T.columns: T2[c]=[('%.2f'%v if abs(v)>=.32 else '  .  ') for v in T[c]]
    T2['h2']=h2.round(2); T2['Ort']=dfc[C22].mean().round(2).values
    print("\n=== 22 madde, k=%d ==="%k); print(T2.to_string())
    print("Phi:",np.round(Phi[np.triu_indices(k,1)],2))
# tek boyutlu cozumden h2<.30 maddeleri de atarak 20 maddelik varyant
L1,P1=efa(Rs,1,'oblimin','minres'); h1=communalities(L1,P1)
low=[C22[i] for i in range(22) if h1[i]<.30]
print("\n22 maddede tek faktorde h2<.30:",low)
C20=[c for c in C22 if c not in low]
print("\n"+"="*104)
print("ADAY MODELLERIN KARSILASTIRMASI (robust WLSMV tipi uyum)")
print("="*104)
rows=[]
def ord_alpha(R): p=R.shape[0]; return p/(p-1)*(1-p/R.sum())
for nm,cs,k in [('31 madde, 4 faktor',cols31,4),
                ('25 madde, 3 faktor',[c for c in cols31 if c not in ['M3','M4','M12','M14','M18','M20']],3),
                ('22 madde, 3 faktor',C22,3),('22 madde, 2 faktor',C22,2),
                ('22 madde, TEK faktor',C22,1),('20 madde, TEK faktor',C20,1)]:
    R,_=polychoric_matrix(dfc[cs]); Rss=smooth(R,0.02); G=gamma_for(cs)
    L,Phi=efa(Rss,k,'oblimin','minres'); b=robust_fit(R,L,Phi,G,N,k)
    h2=communalities(L,Phi); A=np.abs(L); o=np.argsort(-A,1); n_=len(cs); ix=np.arange(n_)
    prim=A[ix,o[:,0]]; sec=A[ix,o[:,1]] if k>1 else np.zeros(n_)
    rows.append(dict(Model=nm, madde=len(cs), k=k, RMSEA=round(b['RMSEA'],3), CFI=round(b['CFI'],3),
        TLI=round(b['TLI'],3), SRMR=round(b['SRMR'],3), alfa=round(ord_alpha(R),3),
        min_yuk=round(prim.min(),2), h2_alti=int((h2<.30).sum()),
        binisik=int(((sec>=.32)&((prim-sec)<.20)).sum()),
        Phi_max=round(np.abs(Phi[np.triu_indices(k,1)]).max(),2) if k>1 else np.nan))
print(pd.DataFrame(rows).to_string(index=False))
pickle.dump({'C22':C22,'C20':C20},open('sets.pkl','wb'))
