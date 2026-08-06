import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
from efa_core import smooth, efa, communalities
d=pickle.load(open('Rc.pkl','rb')); dfc=d['dfc']
def schmid_leiman(R,k):
    """Schmid-Leiman donusumu -> genel + grup faktorleri"""
    L,Phi=efa(smooth(R,0.02),k,'oblimin','minres')
    # 2. duzey: Phi uzerinde tek faktor
    Lh,_=efa(smooth(Phi,0.02),1,'oblimin','minres')
    Lh=Lh.ravel(); Lh=np.clip(Lh,-0.99,0.99)
    U=np.sqrt(np.clip(1-Lh**2,1e-6,None))
    g=L@Lh                      # genel faktor yukleri
    Gr=L*U                      # grup faktor yukleri
    return g,Gr
def indices(g,Gr,R):
    p=len(g)
    cg=g.sum()**2
    cs=(Gr.sum(0)**2).sum()
    h2=g**2+(Gr**2).sum(1)
    unique=(1-h2).sum()
    omega_t=(cg+cs)/(cg+cs+unique)
    omega_h=cg/(cg+cs+unique)
    ECV=(g**2).sum()/((g**2).sum()+(Gr**2).sum())
    return omega_h,omega_t,ECV,omega_h/omega_t
def ord_alpha(R):
    p=R.shape[0]; return p/(p-1)*(1-p/R.sum())
for name,cs in [('25 madde',[c for c in dfc.columns if c not in ['M3','M4','M12','M14','M18','M20']]),
                ('22 madde',[c for c in dfc.columns if c not in ['M3','M4','M12','M14','M18','M20','M11','M17','M29']])]:
    R,_=polychoric_matrix(dfc[cs])
    print("\n=== %s ==="%name)
    print("Sirali (ordinal) alfa = %.3f"%ord_alpha(R))
    for k in [2,3,4]:
        g,Gr=schmid_leiman(R,k)
        oh,ot,ecv,ratio=indices(g,Gr,R)
        print("  bifaktor(k=%d): omega_H=%.3f  omega_total=%.3f  ECV=%.3f  wH/wt=%.3f  genel yuk ort=%.2f"%(
            k,oh,ot,ecv,ratio,np.abs(g).mean()))
    print("  Olcut (Rodriguez, Reise & Haviland 2016): ECV>.70 ve omega_H>.70 -> tek boyutlu ele alinabilir")
