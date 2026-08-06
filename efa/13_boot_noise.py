import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
d=pickle.load(open('Rc.pkl','rb')); df=d['dfc']; cols=d['cols']
rng=np.random.default_rng(5); B=150; n=len(df)
iu=np.triu_indices(31,1)
Rp=[]
for b in range(B):
    idx=rng.integers(0,n,n)
    R,_=polychoric_matrix(df.iloc[idx].reset_index(drop=True))
    Rp.append(R[iu])
Rp=np.array(Rp)
sd=Rp.std(0)
print("Polikorik korelasyonlarin bootstrap SE'si (B=%d):"%B)
print("  ortalama SE=%.4f  medyan=%.4f  90.persentil=%.4f  maks=%.4f"%(sd.mean(),np.median(sd),np.percentile(sd,90),sd.max()))
# Pearson SE karsilastirmasi
Pp=[]
for b in range(B):
    idx=rng.integers(0,n,n); Pp.append(df.iloc[idx].corr().values[iu])
Pp=np.array(Pp); sp=Pp.std(0)
print("Pearson SE:   ortalama=%.4f  medyan=%.4f"%(sp.mean(),np.median(sp)))
print("Polikorik/Pearson SE orani: %.2f"%(sd.mean()/sp.mean()))
ii,jj=iu
o=np.argsort(-sd)[:10]
print("\nEn gurultulu ciftler:")
for k in o: print("  %s-%s SE=%.3f"%(cols[ii[k]],cols[jj[k]],sd[k]))
pickle.dump({'sd':sd,'iu':iu},open('boot.pkl','wb'))
