import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from efa_core import smooth, efa, fit_indices, communalities, map_test
from robust_fit import robust_fit
from final_engine import msa_vec
from polychor import polychoric_matrix
pd.set_option('display.width',220)
g=pickle.load(open('gamma.pkl','rb')); boot=g['boot']; cols31=g['cols']; N=306
iu31=np.triu_indices(31,1)
pair_idx={}
for t,(i,j) in enumerate(zip(*iu31)): pair_idx[(cols31[i],cols31[j])]=t
def gamma_for(cols):
    p=len(cols); sel=[]
    for i in range(p):
        for j in range(i+1,p):
            a,b=cols[i],cols[j]
            sel.append(pair_idx[(a,b)] if (a,b) in pair_idx else pair_idx[(b,a)])
    return len(boot)*0+ np.cov(boot[:,sel].T)*310   # Gamma = n * Cov
d=pickle.load(open('Rc.pkl','rb')); dfc=d['dfc']
SETS={
 '31 madde (tam)': [c for c in cols31],
 '25 madde (tavan etkili 6 madde atildi)': [c for c in cols31 if c not in ['M3','M4','M12','M14','M18','M20']],
 '22 madde (+ iliskisiz M11,M17,M29 atildi)': [c for c in cols31 if c not in ['M3','M4','M12','M14','M18','M20','M11','M17','M29']],
}
for name,cs in SETS.items():
    R,_=polychoric_matrix(dfc[cs]); Rs=smooth(R,0.02); G=gamma_for(cs)
    w=np.linalg.eigvalsh(Rs)[::-1]
    f2,_=map_test(Rs,8)
    print("\n"+"="*100)
    print("%s | KMO-benzeri min MSA=%.3f | ozdeger 1./2.=%.2f | MAP->k=%d"%(
        name,msa_vec(R).min(),w[0]/w[1],np.argmin(f2)))
    print("="*100)
    rows=[]
    for k in range(1,6):
        L,Phi=efa(Rs,k,'oblimin','minres')
        a=fit_indices(Rs,L,k,N,Phi); b=robust_fit(R,L,Phi,G,N,k)
        h2=communalities(L,Phi); A=np.abs(L); o=np.argsort(-A,1); n_=len(cs); ix=np.arange(n_)
        prim=A[ix,o[:,0]]; sec=A[ix,o[:,1]] if k>1 else np.zeros(n_)
        cross=int((( sec>=.32)&((prim-sec)<.20)).sum())
        vc=pd.Series(o[:,0]).value_counts()
        rows.append(dict(k=k, chi2_mv=round(b['T_mv'],1), df=round(b['df_star'],1),
            RMSEA=round(b['RMSEA'],3), CFI=round(b['CFI'],3), TLI=round(b['TLI'],3), SRMR=round(b['SRMR'],3),
            naif_RMSEA=round(a['RMSEA'],3), naif_CFI=round(a['CFI'],3),
            yuk_alti=int((prim<.40).sum()), h2_alti=int((h2<.30).sum()), binisik=cross,
            min_fak=int(vc.min()), Phi_max=round(np.abs(Phi[np.triu_indices(k,1)]).max(),2) if k>1 else 0))
    print(pd.DataFrame(rows).to_string(index=False))
