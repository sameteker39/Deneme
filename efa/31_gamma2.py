import pickle, numpy as np, pandas as pd, sys, time, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
d=pickle.load(open('D.pkl','rb')); dfc=d['dfc']; cols=d['cols']
n=len(dfc); iu=np.triu_indices(31,1); rng=np.random.default_rng(99); B=1000
out=np.zeros((B,len(iu[0]))); t=time.time()
for b in range(B):
    R,_=polychoric_matrix(dfc.iloc[rng.integers(0,n,n)].reset_index(drop=True))
    out[b]=R[iu]
    if (b+1)%200==0: print("  %d/%d (%.0f sn)"%(b+1,B,time.time()-t),flush=True)
pickle.dump({'boot':out,'iu':iu,'cols':cols,'n':n},open('G.pkl','wb'))
print("bitti")
