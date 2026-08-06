import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from items import TXT, ITEM_NO
from polychor import polychoric_matrix
from efa_core import smooth, efa, fit_indices, communalities
from robust_fit import robust_fit
from final_engine import msa_vec
pd.set_option('display.width',250)
g=pickle.load(open('G20.pkl','rb')); boot=g['boot']; c20=g['cols']; N=306
pair={}
for t,(i,j) in enumerate(zip(*g['iu'])): pair[(c20[i],c20[j])]=t
def gamma_for(cs):
    sel=[pair[(cs[i],cs[j])] if (cs[i],cs[j]) in pair else pair[(cs[j],cs[i])]
         for i in range(len(cs)) for j in range(i+1,len(cs))]
    return np.cov(boot[:,sel].T)*310
d=pickle.load(open('K17.pkl','rb')); dfc=d['dfc']; K17=d['keep']; sub=d['sub']
def ord_alpha(R): p=R.shape[0]; return p/(p-1)*(1-p/R.sum())
def sl_indices(R,k):
    L,Phi=efa(smooth(R,0.02),k,'oblimin','minres')
    Lh,_=efa(smooth(Phi,0.02),1,'oblimin','minres'); Lh=np.clip(Lh.ravel(),-.99,.99)
    gen=L@Lh; Gr=L*np.sqrt(np.clip(1-Lh**2,1e-6,None))
    h2=gen**2+(Gr**2).sum(1); cg=gen.sum()**2; cs=(Gr.sum(0)**2).sum(); u=(1-h2).sum()
    return cg/(cg+cs+u), (cg+cs)/(cg+cs+u), (gen**2).sum()/((gen**2).sum()+(Gr**2).sum())
print("="*112); print("ADAY MODELLER (17 madde havuzu, robust WLSMV tipi uyum)"); print("="*112)
rows=[]
CAND=[('17 madde, tek faktor',K17,1),('17 madde, 2 faktor',K17,2),('17 madde, 3 faktor',K17,3),
      ('16 madde (-M16), tek faktor',[c for c in K17 if c!='M16'],1),
      ('16 madde (-M16), 2 faktor',[c for c in K17 if c!='M16'],2),
      ('15 madde (-M16,-M21), 2 faktor',[c for c in K17 if c not in ('M16','M21')],2)]
for nm,cs,k in CAND:
    R,_=polychoric_matrix(dfc[cs]); Rs=smooth(R,0.02); G=gamma_for(cs)
    L,Phi=efa(Rs,k,'oblimin','minres'); b=robust_fit(R,L,Phi,G,N,k)
    h2=communalities(L,Phi); A=np.abs(L); o=np.argsort(-A,1); ix=np.arange(len(cs))
    pri=A[ix,o[:,0]]; sec=A[ix,o[:,1]] if k>1 else np.zeros(len(cs))
    vc=pd.Series(o[:,0]).value_counts()
    rows.append(dict(Model=nm,madde=len(cs),k=k,RMSEA=round(b['RMSEA'],3),CFI=round(b['CFI'],3),
      TLI=round(b['TLI'],3),SRMR=round(b['SRMR'],3),alfa=round(ord_alpha(R),3),
      min_yuk=round(pri.min(),2),h2_alti=int((h2<.30).sum()),
      binisik=int(((sec>=.32)&((pri-sec)<.20)).sum()),
      en_kucuk_fak=int(vc.min()),Phi_max=round(np.abs(Phi[np.triu_indices(k,1)]).max(),2) if k>1 else np.nan))
print(pd.DataFrame(rows).to_string(index=False))
R17,_=polychoric_matrix(dfc[K17])
print("\nBifaktor (Schmid-Leiman) - ozsel tek boyutluluk:")
for k in [2,3]:
    oh,ot,ecv=sl_indices(R17,k)
    print("  k=%d: omega_H=%.3f  omega_total=%.3f  ECV=%.3f   (olcut: ECV>.70 ve omega_H>.70)"%(k,oh,ot,ecv))
