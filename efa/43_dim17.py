import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from items import TXT, ITEM_NO
from polychor import polychoric_matrix
from efa_core import smooth, efa, communalities, map_test
from final_engine import msa_vec
from factor_analyzer.factor_analyzer import calculate_kmo
pd.set_option('display.width',250)
d=pickle.load(open('K20.pkl','rb')); dfc=d['dfc']; sub=d['sub']
keep=[c for c in d['keep'] if c not in ['M11','M17','M29']]
R,_=polychoric_matrix(dfc[keep]); Rs=smooth(R,0.02); p=len(keep)
_,kmo=calculate_kmo(dfc[keep].dropna())
print("17 madde | KMO=%.3f | en dusuk MSA=%.3f | min ozdeger=%.4f"%(kmo,msa_vec(R).min(),np.linalg.eigvalsh(R).min()))
w=np.linalg.eigvalsh(Rs)[::-1]
print("Ozdegerler:",np.round(w[:7],3)," Kaiser>1:",(w>1).sum()," 1./2.=%.2f"%(w[0]/w[1]))
f2,f4=map_test(Rs,7); print("MAP -> k=%d / %d"%(np.argmin(f2),np.argmin(f4)))
lam=w.copy(); ek=[]
for j in range(1,8):
    ref=max(((1+np.sqrt(p/306))**2)*(p-np.sum(lam[:j-1]))/(p-j+1),1); ek.append(lam[j-1]>ref)
print("EKC -> k=%d"%(ek.index(False) if False in ek else 7))
for k in [2,3]:
    L,Phi=efa(Rs,k,'oblimin','minres'); h2=communalities(L,Phi); A=np.abs(L); pri=A.argmax(1)
    print("\n"+"="*110); print("k=%d"%k); print("="*110)
    for f in range(k):
        idx=sorted([i for i in range(p) if pri[i]==f],key=lambda i:-A[i,f])
        print("\n--- F%d (%d madde) ---"%(f+1,len(idx)))
        for i in idx:
            oth="".join("  F%d=%+.2f"%(j+1,L[i,j]) for j in range(k) if j!=f and abs(L[i,j])>=.30)
            print("   %+.2f  %-4s %-55s h2=%.2f%s"%(L[i,f],keep[i],TXT[ITEM_NO[keep[i]]][:55],h2[i],oth))
    print("\nPhi:",np.round(Phi[np.triu_indices(k,1)],2))
pickle.dump({'keep':keep,'R':R,'dfc':dfc[keep],'sub':sub[keep]},open('K17.pkl','wb'))
