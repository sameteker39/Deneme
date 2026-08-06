import numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
from efa_core import smooth, efa, fit_indices, communalities
from robust_fit import robust_fit
rng=np.random.default_rng(321)
p,k,N=15,3,310
L=np.zeros((p,k))
for j in range(k): L[j*5:(j+1)*5,j]=[.80,.75,.70,.65,.60]
Phi=np.array([[1,.5,.4],[.5,1,.45],[.4,.45,1]])
S=L@Phi@L.T; np.fill_diagonal(S,1.0)
cuts=[-2.3,-1.7,-1.1,-0.4]      # agir tavan etkisi (gercek veriye benzer)
X=rng.multivariate_normal(np.zeros(p),S,N)
O=np.digitize(X,cuts)+1.0
df=pd.DataFrame(O,columns=['V%d'%i for i in range(p)])
print("uretilen maddelerin carpikligi: %.2f .. %.2f"%(df.skew().min(),df.skew().max()))
R,_=polychoric_matrix(df)
print("polikorik vs gercek S: ort|fark|=%.3f"%np.abs(R-S)[np.triu_indices(p,1)].mean())
# bootstrap Gamma
B=800; iu=np.triu_indices(p,1); bo=np.zeros((B,len(iu[0])))
for b in range(B):
    idx=rng.integers(0,N,N)
    Rb,_=polychoric_matrix(df.iloc[idx].reset_index(drop=True)); bo[b]=Rb[iu]
G=N*np.cov(bo.T)
print("\n GERCEK MODEL k=3.  Beklenti: k=3'te RMSEA~0, CFI~1\n")
Rs=smooth(R,0.02)
print("%-4s %-34s %s"%("k","NAIF ML (polikorik)","ROBUST (WLSMV tipi)"))
for kk in [1,2,3,4]:
    Lh,Ph=efa(Rs,kk,'oblimin','minres')
    a=fit_indices(Rs,Lh,kk,N,Ph); b=robust_fit(R,Lh,Ph,G,N,kk)
    print("k=%d  RMSEA=%.3f CFI=%.3f SRMR=%.3f   |  RMSEA=%.3f CFI=%.3f TLI=%.3f  (olcek c=%.2f)"%(
        kk,a['RMSEA'],a['CFI'],a['SRMR'],b['RMSEA'],b['CFI'],b['TLI'],b['scale_c']))
