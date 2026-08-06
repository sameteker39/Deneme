"""Sirali veri icin robust (WLSMV tipi) uyum indeksleri.
Muthen (1984), Satorra-Bentler (1994), Asparouhov & Muthen (2010) izleginde:
 - Kestirim: ULS (polikorik korelasyon uzerinde)
 - T = N * (r-sigma)'(r-sigma)
 - U = I - Delta (Delta'Delta)^- Delta'
 - Ortalama-ve-varyans duzeltmeli: df* = tr(UG)^2/tr((UG)^2),  T* = T*df*/tr(UG)
"""
import numpy as np
from scipy.stats import chi2 as _c2

def delta_jacobian(L):
    """d sigma_ij / d lambda_ab  (i<j),  sigma_ij = sum_m l_im l_jm."""
    p,k = L.shape
    iu = np.triu_indices(p,1)
    q = len(iu[0])
    D = np.zeros((q, p*k))
    for t,(i,j) in enumerate(zip(*iu)):
        for b in range(k):
            D[t, i*k+b] += L[j,b]
            D[t, j*k+b] += L[i,b]
    return D

def robust_fit(R, L, Phi, Gamma, N, k):
    """R: gozlenen polikorik; L,Phi: dondurulmus cozum; Gamma: sqrt(N)*r'nin asimptotik kov."""
    p = R.shape[0]; iu = np.triu_indices(p,1)
    C = L@Phi@L.T
    r = R[iu]; s = C[iu]
    e = r - s
    T = N * (e@e)
    # dondurmeden bagimsiz yukleme (Cholesky benzeri) -> Delta
    w,V = np.linalg.eigh(C)
    w = np.clip(w,0,None); o=np.argsort(-w)[:k]
    A = V[:,o]*np.sqrt(w[o])
    D = delta_jacobian(A)
    Dp = np.linalg.pinv(D.T@D, rcond=1e-10)
    U = np.eye(len(r)) - D@Dp@D.T
    UG = U@Gamma
    t1 = np.trace(UG); t2 = np.trace(UG@UG)
    df_nom = ((p-k)**2 - p - k)/2
    df_star = t1*t1/t2 if t2>0 else df_nom
    T_sc  = T*df_nom/t1 if t1>0 else np.nan            # ortalama duzeltmeli (WLSM)
    T_mv  = T*df_star/t1 if t1>0 else np.nan           # ortalama+varyans (WLSMV)
    # bagimsizlik modeli (baseline)
    Db = np.zeros((len(r),0))
    Ub = np.eye(len(r))
    UGb = Ub@Gamma
    t1b = np.trace(UGb); t2b=np.trace(UGb@UGb)
    dfb_star = t1b*t1b/t2b
    Tb = N*(r@r)
    Tb_mv = Tb*dfb_star/t1b
    d_  = max(T_mv-df_star,0); db_ = max(Tb_mv-dfb_star,0)
    rmsea = np.sqrt(d_/(df_star*N)) if df_star>0 else np.nan
    cfi = 1 - d_/max(db_,d_,1e-12)
    tli = ((Tb_mv/dfb_star)-(T_mv/df_star))/((Tb_mv/dfb_star)-1)
    srmr = np.sqrt((e**2).mean())
    return dict(T=T, T_scaled=T_sc, T_mv=T_mv, df=df_nom, df_star=df_star,
                p_mv=1-_c2.cdf(T_mv, max(df_star,1)),
                RMSEA=rmsea, CFI=cfi, TLI=tli, SRMR=srmr, scale_c=t1/df_nom)
