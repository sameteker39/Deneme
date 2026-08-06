import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from purify import purify, diagnose, CRIT
from efa_core import smooth, communalities
d=pickle.load(open('Rc.pkl','rb')); R=d['Rraw']; cols=d['cols']; N=306
res={}
for k in [3,4,5,6]:
    print("\n"+"="*78); print("k = %d  icin yinelemeli madde ayiklama"%k); print("="*78)
    Rf, cf, dd, L, Phi, fi, log = purify(R, cols, k, N)
    res[k]=dict(R=Rf,cols=cf,d=dd,L=L,Phi=Phi,fi=fi,log=log)
    print("  -> kalan %d madde | RMSEA=%.3f CFI=%.3f TLI=%.3f SRMR=%.3f"%(
        len(cf),fi['RMSEA'],fi['CFI'],fi['TLI'],fi['SRMR']))
    print("  faktor buyuklukleri:", dd['Faktor'].value_counts().sort_index().to_dict())
    print("  final min ozdeger (ham polikorik): %.4f"%np.linalg.eigvalsh(Rf).min())
pickle.dump(res, open('purified.pkl','wb'))
