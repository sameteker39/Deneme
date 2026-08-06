import pickle, numpy as np, pandas as pd, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix
from efa_core import smooth
d=pickle.load(open('R.pkl','rb')); df=d['df'].copy(); cols=d['cols']

def collapse(s, minn=5):
    """Alttan yukari dogru seyrek kategorileri birlestir (tavan etkisi yonunde)."""
    v=s.copy()
    while True:
        cnt=v.value_counts().sort_index()
        if len(cnt)<=2 or cnt.min()>=minn: break
        lvl=cnt.idxmin()
        order=sorted(cnt.index)
        i=order.index(lvl)
        tgt=order[i+1] if i+1<len(order) else order[i-1]
        v=v.replace(lvl,tgt)
    return v

dfc=df.copy(); info=[]
for c in cols:
    before=df[c].nunique()
    dfc[c]=collapse(df[c])
    info.append((c,before,dfc[c].nunique(), int(dfc[c].value_counts().min())))
inf=pd.DataFrame(info,columns=['Madde','Once_k','Sonra_k','Min_hucre'])
print("Kategori birlestirilen maddeler:")
print(inf[inf.Once_k!=inf.Sonra_k].to_string(index=False))
print("\nBirlestirilmeyen:", (inf.Once_k==inf.Sonra_k).sum(),"madde")

Rc,_=polychoric_matrix(dfc)
print("\nBirlestirilmis polikorik: min_ozdeger=%.4f cond=%.1f det=%.3e"%(
    np.linalg.eigvalsh(Rc).min(), np.linalg.cond(smooth(Rc)), np.linalg.det(smooth(Rc))))
iu=np.triu_indices(31,1)
print("max |r| = %.3f"%np.abs(Rc[iu]).max())
pickle.dump({'R':smooth(Rc,0.05),'Rraw':Rc,'cols':cols,'dfc':dfc,'df':df},open('Rc.pkl','wb'))
