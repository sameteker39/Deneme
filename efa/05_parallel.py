import pickle, numpy as np, pandas as pd, sys, time
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
from efa_core import smooth

d=pickle.load(open('R.pkl','rb')); df=d['df'].copy(); cols=d['cols']
R=smooth(d['R']); p=R.shape[0]
rng=np.random.default_rng(2024)
NREP=500

def fa_eig(Rm):
    """indirgenmis korelasyon matrisi (kosegen = SMC) ozdegerleri"""
    Ri=np.linalg.inv(smooth(Rm,0.02))
    smc=1-1/np.diag(Ri)
    Rr=Rm.copy(); np.fill_diagonal(Rr,smc)
    return np.linalg.eigvalsh(Rr)[::-1]

obs_pc=np.linalg.eigvalsh(R)[::-1]
obs_fa=fa_eig(R)

V=df.values.copy()
pc=np.zeros((NREP,p)); fa=np.zeros((NREP,p))
t=time.time()
for r in range(NREP):
    P=V.copy()
    for j in range(p):                      # sutun ici permutasyon -> marjinaller korunur
        col=P[:,j]; ok=~np.isnan(col)
        col[ok]=rng.permutation(col[ok]); P[:,j]=col
    Rr,_=polychoric_matrix(pd.DataFrame(P,columns=cols))
    Rr=smooth(Rr)
    pc[r]=np.linalg.eigvalsh(Rr)[::-1]
    fa[r]=fa_eig(Rr)
    if (r+1)%50==0: print("  %d/%d  (%.0f sn)"%(r+1,NREP,time.time()-t),flush=True)

out=pd.DataFrame({
 'k':range(1,p+1),
 'Gozlenen_PC':obs_pc.round(3),'Rastgele_ort_PC':pc.mean(0).round(3),'Rastgele_95_PC':np.percentile(pc,95,0).round(3),
 'Gozlenen_FA':obs_fa.round(3),'Rastgele_ort_FA':fa.mean(0).round(3),'Rastgele_95_FA':np.percentile(fa,95,0).round(3)})
print("\n--- Paralel Analiz (Horn, 500 permutasyon, polikorik) ---")
print(out.head(12).to_string(index=False))
for nm,o,s in [('PC-ortalama',obs_pc,pc.mean(0)),('PC-95%',obs_pc,np.percentile(pc,95,0)),
               ('FA-ortalama',obs_fa,fa.mean(0)),('FA-95%',obs_fa,np.percentile(fa,95,0))]:
    k=int(np.argmax(o<=s)) if (o<=s).any() else p
    print("%-12s -> tutulacak faktor sayisi = %d"%(nm,k))
pickle.dump({'obs_pc':obs_pc,'obs_fa':obs_fa,'pc':pc,'fa':fa}, open('pa.pkl','wb'))
