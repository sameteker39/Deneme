"""
Hızlı tahmin (rapid guessing) eşik yöntemleri.

NT  : Normative Threshold (Wise & Ma, 2012)  -> min(pct*ort_sure, maxsec)
MRTQ: Mixture Response Time Quantile        -> log-RT üzerinde 2 bileşenli
      normal karışım; hızlı-tahmin bileşeninin (düşük ortalamalı) karışım
      ağırlığı w -> eşik = exp( w-inci kuantil(log-RT) ).  (Eker & Gelbal;
      normalmixEM mantığının bağımsız uygulaması.)
"""
import numpy as np

NT_MINN = 30
WINSOR_SD = 2


def winsor_top(x, k=WINSOR_SD):
    x = np.asarray(x, float)
    ok = np.isfinite(x) & (x > 0)
    if ok.sum() < 2:
        return x
    m, s = np.nanmean(x[ok]), np.nanstd(x[ok], ddof=1)
    cap = m + k * s
    return np.where(ok & (x > cap), cap, x)


def nt_threshold(x, pct=0.10, maxsec=10.0, fixed=None, min_n=NT_MINN):
    x = np.asarray(x, float)
    x = x[np.isfinite(x) & (x > 0)]
    if x.size < min_n:
        return np.nan
    if fixed is not None:
        return float(fixed)
    return float(min(pct * x.mean(), maxsec))


def _normmix2_em(lx, max_iter=500, tol=1e-6, seed=0):
    """1B iki bileşenli normal karışım EM. lx: log-RT vektörü."""
    lx = np.asarray(lx, float)
    lx = lx[np.isfinite(lx)]
    n = lx.size
    if n < NT_MINN or np.unique(lx).size < 2:
        return None
    q1, q3 = np.quantile(lx, [0.25, 0.75])
    mu = np.array([q1, q3])
    sd = np.array([lx.std() + 1e-3, lx.std() + 1e-3])
    lam = np.array([0.5, 0.5])
    ll_old = -np.inf
    for _ in range(max_iter):
        # E
        d = np.stack([lam[k] * np.exp(-0.5 * ((lx - mu[k]) / sd[k]) ** 2) /
                      (sd[k] * np.sqrt(2 * np.pi)) for k in range(2)], 1)
        tot = d.sum(1) + 1e-300
        resp = d / tot[:, None]
        ll = np.log(tot).sum()
        # M
        Nk = resp.sum(0) + 1e-9
        lam = Nk / n
        mu = (resp * lx[:, None]).sum(0) / Nk
        var = (resp * (lx[:, None] - mu) ** 2).sum(0) / Nk
        sd = np.sqrt(np.clip(var, 1e-6, None))
        if abs(ll - ll_old) < tol:
            break
        ll_old = ll
    return dict(mu=mu, sd=sd, lam=lam)


def mrtq_threshold(x, min_n=NT_MINN):
    """MRTQ: eşik = exp( w-inci kuantil(log-RT) ), w = hızlı-tahmin bileşen ağırlığı."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x) & (x > 0)]
    if x.size < min_n:
        return np.nan
    lx = np.log(x)
    fit = _normmix2_em(lx)
    if fit is None:
        return np.nan
    rg = int(np.argmin(fit["mu"]))          # düşük ortalamalı bileşen = hızlı tahmin
    w = float(fit["lam"][rg])
    if w < 0.03 or w > 0.97:
        return np.nan
    return float(np.exp(np.quantile(lx, w)))


def classify_rg(fat, items, rule):
    """
    fat: (N,J) süre matrisi (saniye, geçersiz -> nan).
    rule: dict -> {'method':'NT','pct':.10,'maxsec':10,'fixed':None,'winsor':True}
                  veya {'method':'MRTQ','winsor':True}
    Döner: (rgv (N,J) 0/1, thr (J,) eşikler).
    """
    fat = np.asarray(fat, float)
    N, J = fat.shape
    rgv = np.zeros((N, J))
    thr = np.full(J, np.nan)
    winsor = rule.get("winsor", True)
    for j in range(J):
        col = fat[:, j]
        base = winsor_top(col) if winsor else col
        if rule["method"] == "NT":
            thr[j] = nt_threshold(base, rule.get("pct", .10),
                                  rule.get("maxsec", 10.0), rule.get("fixed"))
        else:  # MRTQ
            thr[j] = mrtq_threshold(base)
        t = thr[j]
        rgv[:, j] = np.where(np.isfinite(t) & np.isfinite(col) & (col > 0) & (col <= t), 1, 0)
    return rgv, thr
