import pandas as pd, numpy as np
from scipy import stats
from factor_analyzer.factor_analyzer import calculate_kmo, calculate_bartlett_sphericity

pd.set_option('display.width', 200); pd.set_option('display.max_columns', 50)
SRC='/root/.claude/uploads/e45f5564-4a80-51fe-a001-9fa87c7cfa26/6a6ea3ae-FAfA_veri_31madde.xlsx'
df = pd.read_excel(SRC)
print("Boyut:", df.shape)
print("\n--- Eksik veri ---")
mi = df.isna().sum()
print(mi[mi>0] if mi.sum() else "Yok")
print("Toplam eksik hücre:", int(mi.sum()), " (%.2f%%)" % (100*mi.sum()/df.size))
print("Eksiği olan satır sayısı:", int(df.isna().any(axis=1).sum()))

print("\n--- Değer aralıkları ---")
print("min:", df.min().min(), "max:", df.max().max())
print(pd.Series(np.unique(df.values[~pd.isna(df.values)])).to_list())

d = pd.DataFrame({
 'n': df.notna().sum(),
 'Ort': df.mean().round(2),
 'SS': df.std().round(2),
 'Çarpıklık': df.skew().round(2),
 'Basıklık': df.kurt().round(2),
})
# madde-toplam korelasyonu (düzeltilmiş)
X = df.copy()
tot = X.sum(axis=1)
d['MTK'] = [round(X[c].corr(tot - X[c]),2) for c in X.columns]
print("\n--- Betimsel ---")
print(d)
print("\n|Çarpıklık|>2 veya |Basıklık|>7 olan maddeler:",
      d.index[(d['Çarpıklık'].abs()>2)|(d['Basıklık'].abs()>7)].tolist())
print("MTK < .30 olan maddeler:", d.index[d['MTK']<.30].tolist())

# Mahalanobis (çok değişkenli aykırı)
Z = df.fillna(df.median())
S = np.cov(Z.T); Si = np.linalg.pinv(S); mu = Z.mean().values
md = np.array([ (r-mu)@Si@(r-mu) for r in Z.values ])
crit = stats.chi2.ppf(.999, df.shape[1])
print("\nMahalanobis D2 kritik (p<.001, df=%d) = %.2f | aşan vaka sayısı = %d"
      % (df.shape[1], crit, (md>crit).sum()))

Zc = df.dropna()
print("\nListwise n =", len(Zc))
kmo_item, kmo_all = calculate_kmo(Zc)
chi2, p = calculate_bartlett_sphericity(Zc)
print("KMO (genel) = %.3f" % kmo_all)
print("Bartlett: chi2=%.2f, df=%d, p=%.5f" % (chi2, 31*30/2, p))
print("\nMadde bazlı MSA:")
print(pd.Series(kmo_item, index=Zc.columns).round(3).sort_values().to_string())
