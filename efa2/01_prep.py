import pandas as pd, numpy as np, pickle, sys, warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'.')
from polychor import polychoric_matrix
from factor_analyzer.factor_analyzer import calculate_kmo, calculate_bartlett_sphericity
from final_engine import msa_vec
pd.set_option('display.width',220)
SRC='/root/.claude/uploads/e45f5564-4a80-51fe-a001-9fa87c7cfa26/ff3b842a-ppass.xlsx'
df=pd.read_excel(SRC)
for c in df.columns: df[c]=pd.to_numeric(df[c].astype(str).str.strip(), errors='coerce')
print("Temizlik sonrasi: boyut=%s, eksik=%d (%.2f%%), degerler=%s"%(
    df.shape, df.isna().sum().sum(), 100*df.isna().sum().sum()/df.size,
    sorted(pd.unique(df.values[~pd.isna(df.values)]))))
REV=[3,4,10,11,16]; revc=['M%d'%i for i in REV]
tot=df.sum(axis=1)
print("\n=== Ters kodlama oncesi madde-toplam korelasyonlari ===")
mtk0=pd.Series({c: round(df[c].corr(tot-df[c]),3) for c in df.columns})
print(mtk0.to_string())
print("\nTers kodlanacak:",revc)
for c in revc: df[c]=6-df[c]
tot=df.sum(axis=1)
mtk=pd.Series({c: round(df[c].corr(tot-df[c]),3) for c in df.columns})
print("\n=== Ters kodlama sonrasi ===")
cmp=pd.DataFrame({'Once':mtk0,'Sonra':mtk,'Ters?':[c in revc for c in df.columns]})
print(cmp.to_string())
neg=mtk[mtk<0]
print("\nHala negatif MTK:", neg.to_dict() if len(neg) else "yok -> tum maddeler ayni yonde")
print("MTK<.30:", mtk[mtk<.30].to_dict())
d=pd.DataFrame({'n':df.notna().sum(),'Ort':df.mean().round(2),'SS':df.std().round(2),
  'Carpiklik':df.skew().round(2),'Basiklik':df.kurt().round(2),
  'Tepe%':(100*df.apply(lambda s:s.value_counts().max()/s.notna().sum())).round(1),
  'MinHucre':df.apply(lambda s:int(s.value_counts().reindex([1,2,3,4,5]).fillna(0).min()))})
print("\n=== Betimsel ===\n", d.to_string())
_,kmo=calculate_kmo(df.dropna()); chi2,pv=calculate_bartlett_sphericity(df.dropna())
print("\nKMO=%.3f | Bartlett chi2=%.1f df=%d p=%.2e | listwise n=%d"%(kmo,chi2,21*20/2,pv,len(df.dropna())))
def collapse(s,minn=5):
    v=s.copy()
    while True:
        cnt=v.value_counts().sort_index()
        if len(cnt)<=2 or cnt.min()>=minn: break
        lvl=cnt.idxmin(); o=sorted(cnt.index); i=o.index(lvl)
        v=v.replace(lvl, o[i+1] if i+1<len(o) else o[i-1])
    return v
dfc=df.apply(collapse)
ch=[c for c in df.columns if df[c].nunique()!=dfc[c].nunique()]
print("Kategori birlestirilen:",ch if ch else "yok")
R,cols=polychoric_matrix(dfc)
print("Polikorik: min ozdeger=%.4f  max|r|=%.3f  en dusuk MSA=%.3f"%(
    np.linalg.eigvalsh(R).min(), np.abs(R[np.triu_indices(21,1)]).max(), msa_vec(R).min()))
pickle.dump({'df':df,'dfc':dfc,'R':R,'cols':cols},open('D.pkl','wb'))
