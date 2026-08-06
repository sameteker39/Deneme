import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from items import TXT, ITEM_NO
from efa_core import smooth, efa, communalities, map_test
from final_engine import msa_vec
from factor_analyzer.factor_analyzer import calculate_kmo, calculate_bartlett_sphericity
pd.set_option('display.width',250)
d=pickle.load(open('K20.pkl','rb')); R=d['R']; keep=d['keep']; dfc=d['dfc']; p=len(keep)
Rs=smooth(R,0.02)
_,kmo=calculate_kmo(dfc.dropna()); chi2,pb=calculate_bartlett_sphericity(dfc.dropna())
print("KMO=%.3f | en dusuk MSA=%.3f | Bartlett chi2=%.1f df=%d p<.001"%(kmo,msa_vec(R).min(),chi2,p*(p-1)/2))
w=np.linalg.eigvalsh(Rs)[::-1]
print("\nOzdegerler:",np.round(w[:8],3)," Kaiser>1:",(w>1).sum()," 1./2.=%.2f"%(w[0]/w[1]))
f2,f4=map_test(Rs,8); print("MAP -> k=%d / %d"%(np.argmin(f2),np.argmin(f4)))
J,N=p,306; lam=w.copy(); ek=[]
for j in range(1,9):
    ref=max(((1+np.sqrt(J/N))**2)*(J-np.sum(lam[:j-1]))/(J-j+1),1); ek.append(lam[j-1]>ref)
print("Ampirik Kaiser Olcutu -> k=%d"%(ek.index(False) if False in ek else 8))
for k in [2,3,4]:
    L,Phi=efa(Rs,k,'oblimin','minres'); h2=communalities(L,Phi); A=np.abs(L); pri=A.argmax(1)
    print("\n"+"="*112); print("k=%d"%k); print("="*112)
    for f in range(k):
        idx=sorted([i for i in range(p) if pri[i]==f], key=lambda i:-A[i,f])
        print("\n--- F%d (%d madde) ---"%(f+1,len(idx)))
        for i in idx:
            oth="".join("  F%d=%+.2f"%(j+1,L[i,j]) for j in range(k) if j!=f and abs(L[i,j])>=.30)
            print("   %+.2f  %-4s %-56s h2=%.2f%s"%(L[i,f],keep[i],TXT[ITEM_NO[keep[i]]][:56],h2[i],oth))
    print("\nPhi:",np.round(Phi[np.triu_indices(k,1)],2))
