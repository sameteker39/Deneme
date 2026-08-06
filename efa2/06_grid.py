import pickle, numpy as np, pandas as pd, sys, itertools, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
from efa_core import smooth, efa, communalities, fit_indices
from final_engine import purify_at_k, solution
pd.set_option('display.width',260)
d=pickle.load(open('D.pkl','rb')); R=d['R']; cols=d['cols']; N=202
def ord_alpha(Rm): p=Rm.shape[0]; return p/(p-1)*(1-p/Rm.sum())
rows=[]
for k in [2,3,4,5]:
  for h2t in [0.0,0.20,0.25,0.30]:
    for lt in [0.32,0.40,0.45]:
      for cd in [0.15,0.20]:
        crit=dict(h2=h2t,load=lt,cross_sec=0.32,cross_diff=cd,msa=0.50,minitems=3)
        try: r=purify_at_k(R,cols,k,N,crit=crit,maxdrop=12,verbose=False)
        except Exception: continue
        cs=r['cols']
        if len(cs)<8: continue
        dd=r['d']; vc=dd['Faktor'].value_counts()
        if vc.min()<3: continue                      # <3 maddelik faktor -> gecersiz
        Rk=r['R']; Rs=smooth(Rk,0.02)
        L,Phi=efa(Rs,k,'oblimin','minres'); h2=communalities(L,Phi)
        A=np.abs(L); o=np.argsort(-A,1); ix=np.arange(len(cs))
        pri=A[ix,o[:,0]]; sec=A[ix,o[:,1]]
        fi=fit_indices(Rs,L,k,N,Phi)
        rows.append(dict(k=k,h2t=h2t,lt=lt,cd=cd,madde=len(cs),
          fak=str(vc.sort_index().to_dict()),SRMR=round(fi['SRMR'],3),
          alfa=round(ord_alpha(Rk),3),min_yuk=round(pri.min(),2),min_h2=round(h2.min(),2),
          binisik=int(((sec>=.32)&((pri-sec)<.20)).sum()),
          Phi_max=round(np.abs(Phi[np.triu_indices(k,1)]).max(),2),
          atilan=",".join(x[0] for x in r['log'])))
t=pd.DataFrame(rows).drop_duplicates(subset=['k','atilan'])
t=t.sort_values(['k','madde'],ascending=[True,False])
print("Gecerli cozumler (her faktorde >=3 madde, >=8 madde):",len(t))
print(t.to_string(index=False))
pickle.dump(t,open('grid.pkl','wb'))
