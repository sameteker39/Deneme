import pickle, numpy as np, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
from efa_core import smooth, efa
from scipy.optimize import linear_sum_assignment
d=pickle.load(open('D.pkl','rb')); dfc=d['dfc']; allc=d['cols']
SETS={'C (13 madde)':[c for c in allc if c not in ['M12','M17','M14','M16','M5','M15','M13','M6']],
      'B (15 madde)':[c for c in allc if c not in ['M12','M17','M14','M16','M15','M6']],
      'A (16 madde)':[c for c in allc if c not in ['M12','M17','M20','M6','M15']]}
rng=np.random.default_rng(50); REP=150
def tucker(A,B):
    k=A.shape[1]; M=np.zeros((k,k))
    for i in range(k):
        for j in range(k): M[i,j]=abs(A[:,i]@B[:,j])/np.sqrt((A[:,i]@A[:,i])*(B[:,j]@B[:,j]))
    r,c=linear_sum_assignment(-M); return M[r,c]
for nm,cs in SETS.items():
    sub=dfc[cs]; n=len(sub)
    print("\n=== %s : yari orneklem kararliligi (Tucker phi, %d tekrar) ==="%(nm,REP))
    for k in [1,2,3,4]:
        V=[]
        for _ in range(REP):
            idx=rng.permutation(n); a,b=idx[:n//2],idx[n//2:]
            try:
                Ra,_=polychoric_matrix(sub.iloc[a].reset_index(drop=True))
                Rb,_=polychoric_matrix(sub.iloc[b].reset_index(drop=True))
                La,_=efa(smooth(Ra,0.02),k,'oblimin','minres'); Lb,_=efa(smooth(Rb,0.02),k,'oblimin','minres')
                V.append(tucker(La,Lb))
            except Exception: pass
        V=np.array(V)
        print("  k=%d  ort phi=%.3f | faktorler: %s | phi>=.85 orani=%.2f"%(
            k,V.mean(),np.round(np.sort(V.mean(0))[::-1],3),(V>=.85).mean()))
