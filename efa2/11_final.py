import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
from efa_core import smooth, efa, communalities, fit_indices
from robust_fit import robust_fit
from final_engine import msa_vec
from factor_analyzer.factor_analyzer import calculate_kmo, calculate_bartlett_sphericity
pd.set_option('display.width',260)
g=pickle.load(open('G.pkl','rb')); boot=g['boot']; c21=g['cols']; N=202
pair={}
for t,(i,j) in enumerate(zip(*g['iu'])): pair[(c21[i],c21[j])]=t
def gamma_for(cs):
    sel=[pair[(cs[i],cs[j])] if (cs[i],cs[j]) in pair else pair[(cs[j],cs[i])]
         for i in range(len(cs)) for j in range(i+1,len(cs))]
    return np.cov(boot[:,sel].T)*202
d=pickle.load(open('D.pkl','rb')); dfc=d['dfc']; df=d['df']; allc=d['cols']
C13=[c for c in allc if c not in ['M12','M17','M14','M16','M5','M15','M13','M6']]
R,_=polychoric_matrix(dfc[C13]); Rs=smooth(R,0.02); p=len(C13)
L,Phi=efa(Rs,3,'oblimin','minres'); h2=communalities(L,Phi)
fi=robust_fit(R,L,Phi,gamma_for(C13),N,3)
A=np.abs(L); pri=A.argmax(1)
raw=df[C13]; tot=raw.sum(axis=1)
def oa(Rm): q=Rm.shape[0]; return q/(q-1)*(1-q/Rm.sum())
print("="*104); print("FINAL: 13 madde, 3 FAKTOR (polikorik + ULS/MINRES + oblimin)"); print("="*104)
_,kmo=calculate_kmo(dfc[C13].dropna()); chi2,pv=calculate_bartlett_sphericity(dfc[C13].dropna())
print("KMO=%.3f | en dusuk MSA=%.3f | Bartlett chi2=%.1f df=%d p=%.1e"%(kmo,msa_vec(R).min(),chi2,p*(p-1)/2,pv))
w=np.linalg.eigvalsh(Rs)[::-1]
print("Ozdegerler: %s | Aciklanan toplam varyans: %.1f%%"%(np.round(w[:5],3),100*(L*(L@Phi)).sum()/p))
print("\nUYUM (robust WLSMV tipi): chi2=%.1f df=%.1f | RMSEA=%.3f | CFI=%.3f | TLI=%.3f | SRMR=%.3f"%(
    fi['T_mv'],fi['df_star'],fi['RMSEA'],fi['CFI'],fi['TLI'],fi['SRMR']))
names={0:'F1',1:'F2',2:'F3'}
rows=[]
for f in range(3):
    idx=sorted([i for i in range(p) if pri[i]==f],key=lambda i:-A[i,f])
    for i in idx:
        rows.append(dict(Faktor=names[f],Madde=C13[i],Yuk=round(L[i,f],3),
            **{('F%d'%(j+1)):round(L[i,j],2) for j in range(3)},
            h2=round(h2[i],3),MSA=round(msa_vec(R)[i],3),
            Ort=round(raw[C13[i]].mean(),2),SS=round(raw[C13[i]].std(),2),
            MTK=round(raw[C13[i]].corr(tot-raw[C13[i]]),3)))
print("\n",pd.DataFrame(rows).to_string(index=False))
print("\nFaktorler arasi korelasyon:")
print(pd.DataFrame(Phi.round(3),index=['F1','F2','F3'],columns=['F1','F2','F3']).to_string())
print("\nFaktor bazinda guvenirlik:")
for f in range(3):
    idx=[i for i in range(p) if pri[i]==f]
    Rf=R[np.ix_(idx,idx)]; lf=A[idx,f]
    om=(lf.sum()**2)/((lf.sum()**2)+(1-h2[idx]).sum())
    print("  %s (%d madde): sirali alfa=%.3f  omega=%.3f  ort.yuk=%.2f  maddeler: %s"%(
        names[f],len(idx),oa(Rf),om,lf.mean(),", ".join(C13[i] for i in idx)))
print("  TOPLAM (13 madde): sirali alfa=%.3f"%oa(R))
pd.DataFrame(rows).to_csv('final_13madde_3faktor.csv',index=False)
print("\nDIKKAT: F2'nin tum maddeleri (M3,M4,M10,M11) ters kodlanmis maddelerdir.")
