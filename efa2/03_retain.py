import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
from efa_core import smooth, efa, communalities, map_test
pd.set_option('display.width',250)
d=pickle.load(open('D.pkl','rb')); R=d['R']; cols=d['cols']; dfc=d['dfc']; p=21
Rs=smooth(R,0.02); w=np.linalg.eigvalsh(Rs)[::-1]
print("Ozdegerler:",np.round(w[:9],3))
print("Aciklanan %:",np.round(100*w[:6]/p,2)," Kumulatif:",np.round(100*np.cumsum(w[:6])/p,2))
print("Kaiser>1: %d | 1./2. oran=%.2f"%((w>1).sum(), w[0]/w[1]))
f2,f4=map_test(Rs,8); print("MAP -> k=%d (kare) / %d (4.kuvvet)"%(np.argmin(f2),np.argmin(f4)))
lam=w.copy(); ek=[]
for j in range(1,9):
    ref=max(((1+np.sqrt(p/202))**2)*(p-np.sum(lam[:j-1]))/(p-j+1),1); ek.append(lam[j-1]>ref)
print("Ampirik Kaiser Olcutu -> k=%d"%(ek.index(False) if False in ek else 8))
def fa_eig(Rm):
    Ri=np.linalg.inv(smooth(Rm,0.02)); smc=1-1/np.diag(Ri)
    Rr=Rm.copy(); np.fill_diagonal(Rr,smc); return np.linalg.eigvalsh(Rr)[::-1]
rng=np.random.default_rng(11); B=500; V=dfc.values
opc=w; ofa=fa_eig(Rs); pc=np.zeros((B,p)); fa=np.zeros((B,p))
for b in range(B):
    P=V.copy()
    for j in range(p):
        col=P[:,j]; ok=~np.isnan(col); col[ok]=rng.permutation(col[ok]); P[:,j]=col
    Rb,_=polychoric_matrix(pd.DataFrame(P,columns=cols)); Rb=smooth(Rb,0.02)
    pc[b]=np.linalg.eigvalsh(Rb)[::-1]; fa[b]=fa_eig(Rb)
print("\n=== Paralel Analiz (B=%d, permutasyon, polikorik) ==="%B)
print(pd.DataFrame({'k':range(1,7),'Gozl_PC':opc[:6].round(3),'Rast_ort_PC':pc.mean(0)[:6].round(3),
  'Rast95_PC':np.percentile(pc,95,0)[:6].round(3),'Gozl_FA':ofa[:6].round(3),
  'Rast_ort_FA':fa.mean(0)[:6].round(3),'Rast95_FA':np.percentile(fa,95,0)[:6].round(3)}).to_string(index=False))
for lbl,o,s in [('PC-ort',opc,pc.mean(0)),('PC-95%',opc,np.percentile(pc,95,0)),
                ('FA-ort',ofa,fa.mean(0)),('FA-95%',ofa,np.percentile(fa,95,0))]:
    print("   %-8s -> k=%d"%(lbl,int(np.argmax(o<=s)) if (o<=s).any() else p))
