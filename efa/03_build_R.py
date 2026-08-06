import pandas as pd, numpy as np, sys, time, pickle
sys.path.insert(0,'/home/user/Deneme/efa')
from polychor import polychoric_matrix, nearest_pd

SRC='/root/.claude/uploads/e45f5564-4a80-51fe-a001-9fa87c7cfa26/6a6ea3ae-FAfA_veri_31madde.xlsx'
df = pd.read_excel(SRC)
df['M9'] = 6 - df['M9']                       # ters kodlama
print("M9 ters kodlandi. Yeni MTK:", round(df['M9'].corr(df.sum(axis=1)-df['M9']),3))

# ortalama ikili tam gozlem n
cols=df.columns.tolist(); ns=[]
for i in range(31):
    for j in range(i+1,31):
        ns.append(int((df[cols[i]].notna()&df[cols[j]].notna()).sum()))
print("ikili tam n: min=%d  ort=%.1f  max=%d"%(min(ns),np.mean(ns),max(ns)))

t=time.time()
R, cols = polychoric_matrix(df, verbose=True)
print("\npolikorik matris suresi: %.1f sn"%(time.time()-t))
R = nearest_pd(R)
w=np.linalg.eigvalsh(R)
print("En kucuk ozdeger: %.4f  determinant: %.3e"%(w.min(), np.linalg.det(R)))
Rp = df.corr().values
print("Polikorik vs Pearson ortalama |fark|: %.3f  maks: %.3f"%(
      np.abs(R-Rp)[np.triu_indices(31,1)].mean(), np.abs(R-Rp)[np.triu_indices(31,1)].max()))
pickle.dump({'R':R,'cols':cols,'df':df,'n_mean':float(np.mean(ns))}, open('R.pkl','wb'))
print("\nPolikorik ozdegerler (ilk 12):", np.round(np.linalg.eigvalsh(R)[::-1][:12],3))
