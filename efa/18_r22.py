import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
from efa_core import smooth, efa, fit_indices, communalities, map_test
from factor_analyzer.factor_analyzer import calculate_kmo
from final_engine import msa_vec
pd.set_option('display.width',200)
d=pickle.load(open('Rc.pkl','rb')); dfc=d['dfc']
DROP=['M3','M4','M12','M14','M18','M20','M11','M17','M29']
keep=[c for c in dfc.columns if c not in DROP]
sub=dfc[keep]; R,_=polychoric_matrix(sub); p=len(keep)
print("22 madde: min ozdeger=%.4f cond=%.1f max|r|=%.3f"%(
    np.linalg.eigvalsh(R).min(),np.linalg.cond(R),np.abs(R[np.triu_indices(p,1)]).max()))
_,kmo=calculate_kmo(sub.dropna()); print("KMO=%.3f  min MSA=%.3f"%(kmo,msa_vec(R).min()))
Rs=smooth(R,0.02); w=np.linalg.eigvalsh(Rs)[::-1]
print("Ozdegerler:",np.round(w[:8],3)," Kaiser>1:",(w>1).sum(),
      " | 1./2. oran=%.2f"%(w[0]/w[1]))
f2,f4=map_test(Rs,8); print("MAP -> k=%d / %d"%(np.argmin(f2),np.argmin(f4)))
print("\n--- Uyum (22 madde) ---")
for k in range(1,6):
    L,Phi=efa(Rs,k,'oblimin','minres'); fi=fit_indices(Rs,L,k,306,Phi)
    h2=communalities(L,Phi); A=np.abs(L)
    prim=A.max(1)
    print(" k=%d chi2=%7.1f df=%3d RMSEA=%.3f CFI=%.3f TLI=%.3f SRMR=%.3f | yuk<.40:%d h2<.30:%d"%(
        k,fi['chi2'],fi['df'],fi['RMSEA'],fi['CFI'],fi['TLI'],fi['SRMR'],(prim<.40).sum(),(h2<.30).sum()))
pickle.dump({'R':R,'cols':keep,'sub':sub},open('R22.pkl','wb'))
for k in [1,2]:
    L,Phi=efa(Rs,k,'oblimin','minres'); h2=communalities(L,Phi)
    T=pd.DataFrame(L.round(2),index=keep,columns=['F%d'%(i+1) for i in range(k)]); T['h2']=h2.round(2)
    print("\n== k=%d =="%k); print(T.to_string())
    if k>1: print("Phi:",np.round(Phi[np.triu_indices(k,1)],2))
