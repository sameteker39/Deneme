import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
from efa_core import smooth, efa, communalities, fit_indices
from robust_fit import robust_fit
from final_engine import msa_vec
pd.set_option('display.width',260)
g=pickle.load(open('G.pkl','rb')); boot=g['boot']; c21=g['cols']; N=202
pair={}
for t,(i,j) in enumerate(zip(*g['iu'])): pair[(c21[i],c21[j])]=t
def gamma_for(cs):
    sel=[pair[(cs[i],cs[j])] if (cs[i],cs[j]) in pair else pair[(cs[j],cs[i])]
         for i in range(len(cs)) for j in range(i+1,len(cs))]
    return np.cov(boot[:,sel].T)*202
d=pickle.load(open('D.pkl','rb')); dfc=d['dfc']; allc=d['cols']
def oa(R): p=R.shape[0]; return p/(p-1)*(1-p/R.sum())
S={'A':[c for c in allc if c not in ['M12','M17','M20','M6','M15']],
   'B':[c for c in allc if c not in ['M12','M17','M14','M16','M15','M6']],
   'C':[c for c in allc if c not in ['M12','M17','M14','M16','M5','M15','M13','M6']],
   'D':[c for c in allc if c not in ['M5','M12','M17','M9','M14']]}
CAND=[('Tum 21 madde',allc,1),('Tum 21 madde',allc,2),('Tum 21 madde',allc,3),('Tum 21 madde',allc,4),
      ('A: 16 madde',S['A'],1),('A: 16 madde',S['A'],2),('A: 16 madde',S['A'],3),
      ('B: 15 madde',S['B'],1),('B: 15 madde',S['B'],2),('B: 15 madde',S['B'],3),
      ('C: 13 madde',S['C'],1),('C: 13 madde',S['C'],2),('C: 13 madde',S['C'],3),('C: 13 madde',S['C'],4),
      ('D: 16 madde',S['D'],2)]
rows=[]
for nm,cs,k in CAND:
    R,_=polychoric_matrix(dfc[cs]); Rs=smooth(R,0.02); G=gamma_for(cs)
    L,Phi=efa(Rs,k,'oblimin','minres'); b=robust_fit(R,L,Phi,G,N,k); nai=fit_indices(Rs,L,k,N,Phi)
    h2=communalities(L,Phi); A=np.abs(L); o=np.argsort(-A,1); ix=np.arange(len(cs))
    pri=A[ix,o[:,0]]; sec=A[ix,o[:,1]] if k>1 else np.zeros(len(cs))
    vc=pd.Series(o[:,0]).value_counts()
    rows.append(dict(Model=nm,k=k,madde=len(cs),RMSEA=round(b['RMSEA'],3),CFI=round(b['CFI'],3),
      TLI=round(b['TLI'],3),SRMR=round(b['SRMR'],3),alfa=round(oa(R),3),
      min_yuk=round(pri.min(),2),h2_alti=int((h2<.30).sum()),
      binisik=int(((sec>=.32)&((pri-sec)<.20)).sum()),
      en_kucuk=int(vc.min()),Phi_max=round(np.abs(Phi[np.triu_indices(k,1)]).max(),2) if k>1 else np.nan,
      naif_RMSEA=round(nai['RMSEA'],3)))
print(pd.DataFrame(rows).to_string(index=False))
