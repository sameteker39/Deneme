import pandas as pd, numpy as np, sys, pickle, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from items import ITEM_NO, REVERSE, TXT
from polychor import polychoric_matrix
SRC='/root/.claude/uploads/e45f5564-4a80-51fe-a001-9fa87c7cfa26/6a6ea3ae-FAfA_veri_31madde.xlsx'
df=pd.read_excel(SRC)
tot=df.sum(axis=1)
print("=== Ters madde denetimi (ham veri) ===")
print("%-5s %-4s %-8s %-56s %s"%("Sutun","m#","MTK","Icerik","Durum"))
for c in df.columns:
    n=ITEM_NO[c]; r=df[c].corr(tot-df[c])
    if n in REVERSE:
        st="ZATEN ters kodlanmis" if r>0 else ">>> TERS KODLANACAK <<<"
        print("%-5s m%-3d %+.3f   %-56s %s"%(c,n,r,TXT[n][:56],st))
need=[c for c in df.columns if ITEM_NO[c] in REVERSE and df[c].corr(tot-df[c])<0]
print("\nTers kodlanacak sutun(lar):",need)
for c in need: df[c]=6-df[c]
tot=df.sum(axis=1)
print("Sonrasi MTK:",{c:round(df[c].corr(tot-df[c]),3) for c in need})
print("Tum maddeler pozitif MTK:", all(df[c].corr(tot-df[c])>0 for c in df.columns))
def collapse(s,minn=5):
    v=s.copy()
    while True:
        cnt=v.value_counts().sort_index()
        if len(cnt)<=2 or cnt.min()>=minn: break
        lvl=cnt.idxmin(); order=sorted(cnt.index); i=order.index(lvl)
        v=v.replace(lvl, order[i+1] if i+1<len(order) else order[i-1])
    return v
dfc=df.apply(collapse)
R,cols=polychoric_matrix(dfc)
pickle.dump({'df':df,'dfc':dfc,'R':R,'cols':cols},open('D.pkl','wb'))
print("\nPolikorik hazir. min ozdeger=%.4f"%np.linalg.eigvalsh(R).min())
