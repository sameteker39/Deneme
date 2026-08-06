"""n=202'de GERCEKTEN 3 faktorlu bir yapinin yari-orneklem Tucker phi degeri ne olur?
Gozlenen phi bu referansla karsilastirilmali."""
import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
from efa_core import smooth, efa, communalities
from scipy.optimize import linear_sum_assignment
d=pickle.load(open('D.pkl','rb')); dfc=d['dfc']
C13=[c for c in d['cols'] if c not in ['M12','M17','M14','M16','M5','M15','M13','M6']]
R,_=polychoric_matrix(dfc[C13]); L,Phi=efa(smooth(R,0.02),3,'oblimin','minres')
S=L@Phi@L.T; np.fill_diagonal(S,1.0)
# gozlenen esikler
th=[]
for c in C13:
    v=dfc[c].dropna().values; u=np.sort(np.unique(v))
    cum=np.cumsum([np.mean(v<=x) for x in u])[:-1]*0+np.array([np.mean(v<=x) for x in u])[:-1]
    from scipy.stats import norm
    th.append(norm.ppf(np.clip(cum,1e-4,1-1e-4)))
def tucker(A,B):
    k=A.shape[1]; M=np.zeros((k,k))
    for i in range(k):
        for j in range(k): M[i,j]=abs(A[:,i]@B[:,j])/np.sqrt((A[:,i]@A[:,i])*(B[:,j]@B[:,j]))
    r,c=linear_sum_assignment(-M); return M[r,c]
rng=np.random.default_rng(2718); n=202; p=len(C13); REP=100
print("Referans: gercekten 3 faktorlu model, n=%d, gozlenen yukler ve esikler kullanildi\n"%n)
for k in [1,2,3,4]:
    V=[]
    for _ in range(REP):
        X=rng.multivariate_normal(np.zeros(p),S,n)
        O=np.zeros_like(X)
        for j in range(p): O[:,j]=np.digitize(X[:,j],th[j])+1
        df=pd.DataFrame(O,columns=C13)
        idx=rng.permutation(n); a,b=idx[:n//2],idx[n//2:]
        try:
            Ra,_=polychoric_matrix(df.iloc[a].reset_index(drop=True))
            Rb,_=polychoric_matrix(df.iloc[b].reset_index(drop=True))
            La,_=efa(smooth(Ra,0.02),k,'oblimin','minres'); Lb,_=efa(smooth(Rb,0.02),k,'oblimin','minres')
            V.append(tucker(La,Lb))
        except Exception: pass
    V=np.array(V)
    print("  k=%d  REFERANS ort phi=%.3f (faktorler %s)"%(k,V.mean(),np.round(np.sort(V.mean(0))[::-1],3)))
