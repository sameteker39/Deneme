import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
from efa_core import smooth, efa, communalities
d=pickle.load(open('Rc.pkl','rb')); dfc=d['dfc']
s=pickle.load(open('sets.pkl','rb')); C22=s['C22']
sub=dfc[C22]; n=len(sub); p=len(C22)
rng=np.random.default_rng(4242)

def tucker(A,B):
    """her faktor icin en iyi eslesmeli Tucker uyum katsayisi (isaret/sira serbest)"""
    k=A.shape[1]
    M=np.zeros((k,k))
    for i in range(k):
        for j in range(k):
            M[i,j]=abs(A[:,i]@B[:,j])/np.sqrt((A[:,i]@A[:,i])*(B[:,j]@B[:,j]))
    from scipy.optimize import linear_sum_assignment
    r,c=linear_sum_assignment(-M)
    return M[r,c]

REP=100
print("Yari-orneklem yapi kararliligi (Tucker uyum katsayisi, %d tekrar)"%REP)
print("Olcut: phi>=.95 esdeger, .85-.94 benzer, <.85 kararsiz (Lorenzo-Seva & ten Berge, 2006)\n")
out={}
for k in [1,2,3,4]:
    vals=[]
    for r_ in range(REP):
        idx=rng.permutation(n); a,b=idx[:n//2],idx[n//2:]
        try:
            Ra,_=polychoric_matrix(sub.iloc[a].reset_index(drop=True))
            Rb,_=polychoric_matrix(sub.iloc[b].reset_index(drop=True))
            La,_=efa(smooth(Ra,0.02),k,'oblimin','minres')
            Lb,_=efa(smooth(Rb,0.02),k,'oblimin','minres')
            vals.append(tucker(La,Lb))
        except Exception: pass
    V=np.array(vals)
    out[k]=V
    print("k=%d  ortalama phi=%.3f | faktor bazinda ort: %s | phi>=.85 orani=%.2f  | en zayif faktor ort=%.3f"%(
        k, V.mean(), np.round(np.sort(V.mean(0))[::-1],3), (V>=.85).mean(), V.mean(0).min()))
pickle.dump(out,open('stab.pkl','wb'))
