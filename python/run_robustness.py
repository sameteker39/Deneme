"""
DUYARLILIK ANALİZİ — hızlı tahmin eşik kurallarının parametre kestirimine etkisi.

Kurallar:
  NT10 (2-SS budamalı, ANA) | NT10 (budasız) | NT5 (budamalı) | Sabit 5 sn | MRTQ

Her kural için: RG oranı, RG geçerliği (doğruluk|RG vs |çaba), ve üç modelin
(Standart / EM-IRT / IRTree-ACC) madde+yetenek parametreleri; ANA (NT10) ile
korelasyon ve ortalama ayırt edicilik. mirt yerine python/irt_est.py (MML-EM).

Çıktı: outputs_python/T5_duyarlilik.csv
"""
import os
import numpy as np
import pandas as pd
import pyreadr
from irt_est import fit_2pl_mml, fit_irtree_2d
from thresholds import classify_rg

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "outputs_python")
os.makedirs(OUT, exist_ok=True)
GECERLI, MAX_TIME = [1, 2, 3, 4], 1e5

# ---- veri ----
d = list(pyreadr.read_r(os.path.join(ROOT, "data", "CALISMA2.Rdata")).values())[0]
keys = pd.read_csv(os.path.join(ROOT, "data", "recovered_keys.csv"))
items = keys["madde_id"].tolist()
corr = dict(zip(keys["madde_id"], keys["dogru_kod"]))
acc, fat = {}, {}
for it in items:
    r = pd.to_numeric(d[it], errors="coerce")
    acc[it] = np.where(r.isin(GECERLI), (r == corr[it]).astype(float), np.nan)
    t = pd.to_numeric(d[it + "_R"], errors="coerce")
    fat[it] = np.where(np.isfinite(t) & (t > 0) & (t < MAX_TIME), t, np.nan)
acc = pd.DataFrame(acc).to_numpy(float)
fat = pd.DataFrame(fat).to_numpy(float)

# Standart model (eşikten bağımsız) bir kez
std = fit_2pl_mml(acc, n_grid=41, max_iter=300)

RULES = {
    "NT10 (budamalı, ANA)": dict(method="NT", pct=.10, maxsec=10, winsor=True),
    "NT10 (budasız)":       dict(method="NT", pct=.10, maxsec=10, winsor=False),
    "NT5 (budamalı)":       dict(method="NT", pct=.05, maxsec=10, winsor=True),
    "Sabit 5 sn":           dict(method="NT", fixed=5, winsor=False),
    "MRTQ":                 dict(method="MRTQ", winsor=True),
}

def corr(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[ok], y[ok])[0, 1])

main_em = None
rows = []
for name, rule in RULES.items():
    print(f">>> {name} ...")
    rgv, thr = classify_rg(fat, items, rule)
    eng = np.isfinite(fat) & (fat > 0)
    rg_rate = rgv[eng].mean()
    both = eng & np.isfinite(acc)
    acc_rg = acc[both & (rgv == 1)].mean() if (both & (rgv == 1)).any() else np.nan
    acc_sol = acc[both & (rgv == 0)].mean()
    # modeller
    E = acc.copy(); E[rgv == 1] = np.nan
    X1 = np.where(eng, rgv, np.nan)
    X2 = np.where(rgv == 0, acc, np.nan)
    em = fit_2pl_mml(E, n_grid=41, max_iter=300)
    tree = fit_irtree_2d(X1, X2, n_grid=25, max_iter=150)
    if main_em is None:
        main_em = em                      # ANA = ilk kural (NT10 budamalı)
    rows.append(dict(
        kural=name,
        RG_orani=round(rg_rate, 3),
        RG_dogruluk=round(acc_rg, 3),
        caba_dogruluk=round(acc_sol, 3),
        a_ort_EMIRT=round(np.mean(em["a"]), 3),
        b_ort_EMIRT=round(np.mean(em["b"]), 3),
        a_ort_IRTreeACC=round(np.mean(tree["a2"]), 3),
        rho_RG_ACC=round(tree["rho"], 3),
        a_kor_ANA=round(corr(main_em["a"], em["a"]), 3),
        b_kor_ANA=round(corr(main_em["b"], em["b"]), 3),
        theta_kor_ANA=round(corr(main_em["theta"], em["theta"]), 3),
    ))

tab = pd.DataFrame(rows)
tab.to_csv(os.path.join(OUT, "T5_duyarlilik.csv"), index=False)
pd.set_option("display.width", 200)
print("\n=== DUYARLILIK TABLOSU ===")
print(tab.to_string(index=False))
