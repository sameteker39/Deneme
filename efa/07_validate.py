import numpy as np, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from efa_core import efa, fit_indices, smooth
rng=np.random.default_rng(7)
# bilinen 3 faktorlu model, 12 madde, temiz yapi
p,k=12,3
L=np.zeros((p,k))
for j in range(k): L[j*4:(j+1)*4, j]=[.75,.70,.65,.60]
Phi=np.array([[1,.4,.3],[.4,1,.35],[.3,.35,1]])
S=L@Phi@L.T; np.fill_diagonal(S,1.0)
N=500
X=rng.multivariate_normal(np.zeros(p), S, N)
R=np.corrcoef(X.T)
print("dogru model kosullari: k=3 icin CFI~1, RMSEA~0 beklenir\n")
for kk in [1,2,3,4]:
    Lh,Ph=efa(R,kk,'oblimin','minres')
    fi=fit_indices(R,Lh,kk,N,Ph)
    print("k=%d  chi2=%7.1f df=%3d  RMSEA=%.3f  CFI=%.3f  TLI=%.3f  SRMR=%.3f"%(
        kk,fi['chi2'],fi['df'],fi['RMSEA'],fi['CFI'],fi['TLI'],fi['SRMR']))
print("\nR kosul sayisi=%.1f  det=%.4f  min ozdeger=%.4f"%(
    np.linalg.cond(R), np.linalg.det(R), np.linalg.eigvalsh(R).min()))
