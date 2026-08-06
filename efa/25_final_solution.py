import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
from efa_core import smooth, efa, fit_indices, communalities
from robust_fit import robust_fit
from final_engine import msa_vec
from factor_analyzer.factor_analyzer import calculate_kmo, calculate_bartlett_sphericity
pd.set_option('display.width',220)
g=pickle.load(open('gamma.pkl','rb')); boot=g['boot']; cols31=g['cols']; N=306
iu31=np.triu_indices(31,1); pair={}
for t,(i,j) in enumerate(zip(*iu31)): pair[(cols31[i],cols31[j])]=t
def gamma_for(cs):
    sel=[pair[(cs[i],cs[j])] if (cs[i],cs[j]) in pair else pair[(cs[j],cs[i])]
         for i in range(len(cs)) for j in range(i+1,len(cs))]
    return np.cov(boot[:,sel].T)*310
d=pickle.load(open('Rc.pkl','rb')); dfc=d['dfc']; dforig=d['df']
C22=pickle.load(open('sets.pkl','rb'))['C22']
R,_=polychoric_matrix(dfc[C22]); Rs=smooth(R,0.02); G=gamma_for(C22)
L,Phi=efa(Rs,1,'oblimin','minres'); h2=communalities(L,Phi); l=L.ravel()
fi=robust_fit(R,L,Phi,G,N,1); nai=fit_indices(Rs,L,1,N,Phi)
raw=dforig[C22]
tot=raw.sum(axis=1)
mtk=[raw[c].corr(tot-raw[c]) for c in C22]
msa=msa_vec(R)
def a_drop(cs):
    Rr,_=polychoric_matrix(dfc[cs]); p=len(cs); return p/(p-1)*(1-p/Rr.sum())
alpha=len(C22)/(len(C22)-1)*(1-len(C22)/R.sum())
om=(l.sum()**2)/((l.sum()**2)+(1-h2).sum())
T=pd.DataFrame({'Madde':C22,'Ort':raw.mean().round(2).values,'SS':raw.std().round(2).values,
 'Carpiklik':raw.skew().round(2).values,'Basiklik':raw.kurt().round(2).values,
 'Faktor_yuku':l.round(3),'h2':h2.round(3),'MSA':msa.round(3),
 'Madde_toplam_r':np.round(mtk,3),
 'Alfa_madde_atilirsa':[round(a_drop([x for x in C22 if x!=c]),3) for c in C22]})
print("="*112); print("FINAL COZUM: 22 madde, TEK FAKTOR (polikorik + ULS/MINRES)"); print("="*112)
print(T.to_string(index=False))
kmo_i,kmo=calculate_kmo(dfc[C22].dropna()); chi2b,pb=calculate_bartlett_sphericity(dfc[C22].dropna())
w=np.linalg.eigvalsh(Rs)[::-1]
print("\nOrneklem uygunlugu: KMO=%.3f | Bartlett chi2=%.1f, df=%d, p<.001"%(kmo,chi2b,22*21/2))
print("Ozdegerler: %s"%np.round(w[:5],3))
print("Aciklanan varyans (tek faktor): %.1f%%  | 1./2. ozdeger orani=%.2f"%(100*(l**2).sum()/22, w[0]/w[1]))
print("\nUYUM IYILIGI (robust, WLSMV tipi):")
print("  chi2(mv)=%.1f, df=%.1f, p=%.3f"%(fi['T_mv'],fi['df_star'],fi['p_mv']))
print("  RMSEA=%.3f | CFI=%.3f | TLI=%.3f | SRMR=%.3f"%(fi['RMSEA'],fi['CFI'],fi['TLI'],fi['SRMR']))
print("  (duzeltilmemis naif ML referansi: RMSEA=%.3f CFI=%.3f -- asiri karamsar)"%(nai['RMSEA'],nai['CFI']))
print("\nGUVENIRLIK: sirali alfa=%.3f | McDonald omega=%.3f"%(alpha,om))
print("Yuk araligi: %.2f - %.2f | tum h2>=%.2f"%(l.min(),l.max(),h2.min()))
DROP=[c for c in cols31 if c not in C22]
print("\nATILAN 9 MADDE:",", ".join(DROP))
T.to_csv('final_22_madde.csv',index=False)
