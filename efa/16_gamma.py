"""Polikorik korelasyonlarin asimptotik kovaryans matrisi (Gamma) - bootstrap."""
import pickle, numpy as np, pandas as pd, sys, time, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
d=pickle.load(open('Rc.pkl','rb')); dfc=d['dfc']; cols=d['cols']
n=len(dfc); p=len(cols); iu=np.triu_indices(p,1)
rng=np.random.default_rng(99); B=1000
out=np.zeros((B,len(iu[0])))
t=time.time()
for b in range(B):
    idx=rng.integers(0,n,n)
    R,_=polychoric_matrix(dfc.iloc[idx].reset_index(drop=True))
    out[b]=R[iu]
    if (b+1)%100==0: print("  %d/%d (%.0f sn)"%(b+1,B,time.time()-t),flush=True)
pickle.dump({'boot':out,'iu':iu,'cols':cols,'n':n},open('gamma.pkl','wb'))
print("bitti. Gamma boyutu:",out.shape)
