import pickle, numpy as np, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
from efa_core import smooth, efa
from scipy.optimize import linear_sum_assignment
d=pickle.load(open('K17.pkl','rb')); dfc=d['dfc'][[c for c in d['keep'] if c not in('M16','M21')]]
rng=np.random.default_rng(999); n=len(dfc); REP=150
def tucker(A,B):
    k=A.shape[1]; M=np.zeros((k,k))
    for i in range(k):
        for j in range(k): M[i,j]=abs(A[:,i]@B[:,j])/np.sqrt((A[:,i]@A[:,i])*(B[:,j]@B[:,j]))
    r,c=linear_sum_assignment(-M); return M[r,c]
print("15 MADDE (-M16,-M21) kararlilik (Tucker phi, %d tekrar)\n"%REP)
for k in [1,2,3]:
    V=[]
    for _ in range(REP):
        idx=rng.permutation(n); a,b=idx[:n//2],idx[n//2:]
        try:
            Ra,_=polychoric_matrix(dfc.iloc[a].reset_index(drop=True))
            Rb,_=polychoric_matrix(dfc.iloc[b].reset_index(drop=True))
            La,_=efa(smooth(Ra,0.02),k,'oblimin','minres'); Lb,_=efa(smooth(Rb,0.02),k,'oblimin','minres')
            V.append(tucker(La,Lb))
        except Exception: pass
    V=np.array(V)
    print("k=%d  ort phi=%.3f | faktorler: %s | phi>=.85 orani=%.2f"%(
        k,V.mean(),np.round(np.sort(V.mean(0))[::-1],3),(V>=.85).mean()))
