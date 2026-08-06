import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
from efa_core import smooth, efa, fit_indices, communalities, map_test
from factor_analyzer.factor_analyzer import calculate_kmo
pd.set_option('display.width',200)
d=pickle.load(open('Rc.pkl','rb')); dfc=d['dfc']; cols=d['cols']
DROP=['M3','M4','M12','M14','M18','M20']
keep=[c for c in cols if c not in DROP]
sub=dfc[keep]
R,_=polychoric_matrix(sub)
print("25 madde polikorik: min ozdeger=%.4f  cond=%.1f  det=%.3e"%(
    np.linalg.eigvalsh(R).min(), np.linalg.cond(R), np.linalg.det(R)))
print("max |r| = %.3f"%np.abs(R[np.triu_indices(25,1)]).max())
kmo_i,kmo=calculate_kmo(sub.dropna()); print("KMO = %.3f"%kmo)
w=np.linalg.eigvalsh(smooth(R,0.02))[::-1]
print("\nOzdegerler:", np.round(w[:10],3), " Kaiser>1:",(w>1).sum())
f2,f4=map_test(smooth(R,0.02),10)
print("MAP(kare) min -> k=%d | MAP(4.kuvvet) -> k=%d"%(np.argmin(f2),np.argmin(f4)))
print("\n--- Uyum (25 madde, polikorik ULS) ---")
for k in range(1,7):
    L,Phi=efa(smooth(R,0.02),k,'oblimin','minres'); fi=fit_indices(smooth(R,0.02),L,k,306,Phi)
    print(" k=%d chi2=%7.1f df=%3d RMSEA=%.3f CFI=%.3f TLI=%.3f SRMR=%.3f"%(
        k,fi['chi2'],fi['df'],fi['RMSEA'],fi['CFI'],fi['TLI'],fi['SRMR']))
pickle.dump({'R':R,'cols':keep,'sub':sub,'dfc':dfc},open('R25.pkl','wb'))
