import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
from efa_core import smooth, efa, communalities
from final_engine import msa_vec
pd.set_option('display.width',250)
d=pickle.load(open('D.pkl','rb')); dfc=d['dfc']; allc=d['cols']
CAND={
 'A: k=3, 16 madde': ([c for c in allc if c not in ['M12','M17','M20','M6','M15']],3),
 'B: k=3, 15 madde': ([c for c in allc if c not in ['M12','M17','M14','M16','M15','M6']],3),
 'C: k=3, 13 madde': ([c for c in allc if c not in ['M12','M17','M14','M16','M5','M15','M13','M6']],3),
 'D: k=2, 16 madde': ([c for c in allc if c not in ['M5','M12','M17','M9','M14']],2),
}
for nm,(cs,k) in CAND.items():
    R,_=polychoric_matrix(dfc[cs]); Rs=smooth(R,0.02); p=len(cs)
    L,Phi=efa(Rs,k,'oblimin','minres'); h2=communalities(L,Phi); A=np.abs(L); pri=A.argmax(1)
    print("\n"+"="*96); print("%s  (KMO-min MSA=%.3f)"%(nm,msa_vec(R).min())); print("="*96)
    for f in range(k):
        idx=sorted([i for i in range(p) if pri[i]==f],key=lambda i:-A[i,f])
        print(" --- F%d (%d madde) ---"%(f+1,len(idx)))
        for i in idx:
            oth="".join("  F%d=%+.2f"%(j+1,L[i,j]) for j in range(k) if j!=f and abs(L[i,j])>=.30)
            print("    %+.2f  %-4s  h2=%.2f%s"%(L[i,f],cs[i],h2[i],oth))
    print(" Phi:",np.round(Phi[np.triu_indices(k,1)],2))
