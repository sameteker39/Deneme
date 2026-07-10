"""
ÖN SONUÇLAR — IRTree vs EM-IRT vs Standart 2PL (Python MML-EM, mirt olmadan).

Amaç: mirt kurulamayan ortamlarda üç modelin madde (a,b) ve yetenek (theta)
parametrelerini kestirip karşılaştırmalı ön sonuç üretmek. Nihai analiz için
R/calisma3_irtree_vs_emirt.R (mirt) kullanılmalıdır.

Kestiriciler python/irt_est.py içinde; ikisi de simülasyonla doğrulanmıştır
(2PL: a kor≈0.99, b kor≈1.00; IRTree: rho ve a/b geri kazanımı ✓).

Girdi:  data/CALISMA2.Rdata , data/recovered_keys.csv
Çıktı:  outputs_python/T1_madde_parametreleri.csv
        outputs_python/T3_madde_karsilastirma.csv
        outputs_python/T4_theta_karsilastirma.csv

Gerekli: pyreadr, pandas, numpy   (opsiyonel: matplotlib)
"""
import os
import numpy as np
import pandas as pd
import pyreadr
from irt_est import fit_2pl_mml, fit_irtree_2d

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data", "CALISMA2.Rdata")
KEYS = os.path.join(ROOT, "data", "recovered_keys.csv")
OUT = os.path.join(ROOT, "outputs_python")
os.makedirs(OUT, exist_ok=True)

GECERLI = [1, 2, 3, 4]
MAX_TIME = 1e5  # gecerli RT ust siniri (sentinel 9999999999 haric)

# ---- veri hazırlama ----
d = list(pyreadr.read_r(DATA).values())[0]
keys = pd.read_csv(KEYS)
items = keys["madde_id"].tolist()
corr = dict(zip(keys["madde_id"], keys["dogru_kod"]))

acc, fat = {}, {}
for it in items:
    r = pd.to_numeric(d[it], errors="coerce")
    acc[it] = np.where(r.isin(GECERLI), (r == corr[it]).astype(float), np.nan)
    t = pd.to_numeric(d[it + "_R"], errors="coerce")
    fat[it] = np.where(np.isfinite(t) & (t > 0) & (t < MAX_TIME), t, np.nan)
acc = pd.DataFrame(acc)
fat = pd.DataFrame(fat)

# ---- NT10 (madde içi 2-SS winsorize) ----
def winsor_top(x, k=2):
    x = x.astype(float)
    ok = np.isfinite(x) & (x > 0)
    m, s = np.nanmean(x[ok]), np.nanstd(x[ok], ddof=1)
    return np.where(ok & (x > m + k * s), m + k * s, x)

rgv = np.zeros(acc.shape)
for j, it in enumerate(items):
    fx = winsor_top(fat[it].values)
    v = fx[np.isfinite(fx) & (fx > 0)]
    thr = min(0.10 * v.mean(), 10.0) if len(v) >= 30 else np.nan
    f = fat[it].values
    rgv[:, j] = np.where(np.isfinite(thr) & np.isfinite(f) & (f > 0) & (f <= thr), 1, 0)

A = acc.to_numpy(float)                 # Standart: tüm puanlı yanıtlar
E = A.copy(); E[rgv == 1] = np.nan       # EM-IRT: hızlı tahmin -> NA
X1 = np.where(np.isfinite(fat.to_numpy()) & (fat.to_numpy() > 0), rgv, np.nan)  # RG düğümü
X2 = np.where(rgv == 0, A, np.nan)       # ACC düğümü (effortful)

# ---- modeller ----
print(">>> Standart 2PL ..."); std = fit_2pl_mml(A, n_grid=41, max_iter=300)
print(">>> EM-IRT ...");       em = fit_2pl_mml(E, n_grid=41, max_iter=300)
print(">>> IRTree 2D ...");    tree = fit_irtree_2d(X1, X2, n_grid=25, max_iter=200)
print(f"    rho(RG,ACC) = {tree['rho']:.3f}")

# ---- tablolar ----
ip = pd.DataFrame({"madde": items,
                   "a_std": std["a"], "b_std": std["b"],
                   "a_em": em["a"], "b_em": em["b"],
                   "a_tree": tree["a2"], "b_tree": tree["b2"]})
ip.round(4).to_csv(os.path.join(OUT, "T1_madde_parametreleri.csv"), index=False)

def cmp(x, y, lab):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    return dict(karsilastirma=lab, n=int(ok.sum()),
                r=round(float(np.corrcoef(x[ok], y[ok])[0, 1]), 3),
                ort_fark=round(float(np.mean(y[ok] - x[ok])), 3),
                RMSD=round(float(np.sqrt(np.mean((y[ok] - x[ok]) ** 2))), 3))

pd.DataFrame([
    cmp(ip.b_std, ip.b_em, "b: Standart->EM-IRT"),
    cmp(ip.b_std, ip.b_tree, "b: Standart->IRTree(ACC)"),
    cmp(ip.b_em, ip.b_tree, "b: EM-IRT->IRTree(ACC)"),
    cmp(ip.a_std, ip.a_em, "a: Standart->EM-IRT"),
    cmp(ip.a_std, ip.a_tree, "a: Standart->IRTree(ACC)"),
    cmp(ip.a_em, ip.a_tree, "a: EM-IRT->IRTree(ACC)"),
]).to_csv(os.path.join(OUT, "T3_madde_karsilastirma.csv"), index=False)

th = pd.DataFrame({"theta_std": std["theta"], "theta_em": em["theta"], "theta_tree": tree["theta2"]})
pd.DataFrame([
    cmp(th.theta_std, th.theta_em, "theta: Standart->EM-IRT"),
    cmp(th.theta_std, th.theta_tree, "theta: Standart->IRTree(ACC)"),
    cmp(th.theta_em, th.theta_tree, "theta: EM-IRT->IRTree(ACC)"),
]).to_csv(os.path.join(OUT, "T4_theta_karsilastirma.csv"), index=False)

print("Tamam. Çıktılar:", OUT)
